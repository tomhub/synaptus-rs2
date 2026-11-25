import argparse
import hashlib
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RS1_ROOT = ROOT.parent / "synaptus-rs"
OUTPUT = Path(__file__).resolve().parent
TOLERANCE = 1.0e-11
LANGUAGES = {
    "rust": ("Rust rs2 PSM", "#0c8177"),
    "rust_v1": ("Rust rs PSM", "#df6848"),
    "octave": ("Octave", "#386cb0"),
    "matlab": ("MATLAB", "#386cb0"),
    "python": ("Python / NumPy", "#b08028"),
}


def command_output(command, cwd=ROOT, env=None):
    return subprocess.run(command, cwd=cwd, env=env, check=True, text=True, capture_output=True).stdout.strip()


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cpu_name():
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "not reported"


def quantile(values, probability):
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def describe(values):
    return {
        "samples_sec": values,
        "mean_sec": statistics.mean(values),
        "median_sec": statistics.median(values),
        "q1_sec": quantile(values, 0.25),
        "q3_sec": quantile(values, 0.75),
        "sd_sec": statistics.stdev(values),
        "min_sec": min(values),
        "max_sec": max(values),
    }


def read_samples(path, expected_count):
    values = [float(value) for value in path.read_text().strip().split(",") if value]
    if len(values) != expected_count or not all(math.isfinite(value) and value > 0 for value in values):
        raise RuntimeError(f"Invalid timing samples in {path}")
    return values


def run_drivers(engine, warmups, repeats, temp_dir, rayon_threads):
    env = os.environ.copy()
    env.update({
        "RAYON_NUM_THREADS": str(rayon_threads),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    })
    binary_name = "speed_benchmark.exe" if os.name == "nt" else "speed_benchmark"
    rust_outputs = {}
    for key, project in (("rust", ROOT), ("rust_v1", RS1_ROOT)):
        output_dir = temp_dir / key
        output_dir.mkdir()
        subprocess.run(["cargo", "build", "--release", "--example", "speed_benchmark"], cwd=project, env=env, check=True)
        binary = project / "target/release/examples" / binary_name
        subprocess.run([str(binary), str(ROOT / "data"), str(output_dir), str(warmups), str(repeats)], cwd=project, env=env, check=True)
        rust_outputs[key] = output_dir

    matlab_script = "psm_benchmark_matlab" if engine == "matlab" else "psm_benchmark"
    output_key = engine
    matlab_path = lambda path: str(path.resolve()).replace("\\", "/").replace("'", "''")
    expression = (
        f"addpath('{matlab_path(OUTPUT)}'); "
        f"{matlab_script}('{matlab_path(ROOT / 'data')}','{matlab_path(temp_dir)}',{warmups},{repeats})"
    )
    if engine == "matlab":
        subprocess.run(["matlab", "-batch", expression], cwd=ROOT, env=env, check=True)
    else:
        subprocess.run(["octave", "--quiet", "--no-gui", "--eval", expression], cwd=ROOT, env=env, check=True)
    subprocess.run([
        sys.executable, "-I", str(OUTPUT / "python_psm_benchmark.py"),
        str(ROOT / "data"), str(temp_dir), str(warmups), str(repeats),
    ], cwd=ROOT, env=env, check=True)

    sample_files = {
        "rust": rust_outputs["rust"] / "rust_samples.csv",
        "rust_v1": rust_outputs["rust_v1"] / "rust_samples.csv",
        output_key: temp_dir / f"{output_key}_samples.csv",
        "python": temp_dir / "python_samples.csv",
    }
    return {key: read_samples(path, repeats) for key, path in sample_files.items()}


def check_images(engine, temp_dir):
    images = {}
    image_paths = {
        "rust": temp_dir / "rust/rust_result.c64",
        "rust_v1": temp_dir / "rust_v1/rust_result.c64",
        engine: temp_dir / f"{engine}_result.c64",
        "python": temp_dir / "python_result.c64",
    }
    for key, path in image_paths.items():
        raw = np.fromfile(path, dtype="<f8")
        if raw.size % 2:
            raise RuntimeError(f"{key} output has a truncated complex value")
        image = raw[0::2] + 1j * raw[1::2]
        if not np.all(np.isfinite(image)):
            raise RuntimeError(f"{key} output contains NaN or infinity")
        images[key] = image
    reference = images["rust"]
    peak = float(np.max(np.abs(reference)))
    if peak == 0.0:
        raise RuntimeError("rs2 Rust output peak is zero")
    agreement = {}
    for key, image in images.items():
        if image.shape != reference.shape:
            raise RuntimeError(f"{key} output shape {image.shape} does not match Rust {reference.shape}")
        error = float(np.max(np.abs(image - reference)) / peak)
        agreement[key] = {"points": int(image.size), "max_abs_error_over_rs2_peak": error}
        if error > TOLERANCE:
            raise RuntimeError(f"{key} output fails the {TOLERANCE:g} normalized-to-rs2-peak tolerance: {error:g}")
    return agreement, int(reference.size)


def collect_environment(engine, rayon_threads):
    if engine == "matlab":
        runtime_version = command_output(["matlab", "-batch", "disp(version)"]).splitlines()[-1]
    else:
        runtime_version = command_output(["octave", "--version"]).splitlines()[0]
    try:
        revision = command_output(["git", "rev-parse", "HEAD"])
    except subprocess.CalledProcessError:
        revision = "unavailable"
    return {
        "os": platform.platform(),
        "cpu": cpu_name(),
        "logical_cpus": os.cpu_count(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        f"{engine}_version": runtime_version,
        "rustc": command_output(["rustc", "--version"]),
        "cargo": command_output(["cargo", "--version"]),
        "parallelism": f"Rust rs and rs2 use Rayon={rayon_threads} workers each; MATLAB/Octave and common NumPy math thread env vars are limited to 1",
        "rayon_threads_per_rust_implementation": rayon_threads,
        "git_revision": revision,
    }


def build_report(report):
    data = json.dumps(report, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    html = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light">
<title>Synaptus rs2 baseline PSM | performance</title>
<style>
:root{color-scheme:light;--ink:#172a2b;--muted:#586a69;--paper:#f7f8f4;--white:#fff;--line:#d5ded8;--green:#0c8177;--orange:#df6848;--blue:#386cb0;--lime:#d9e86a;--shadow:0 12px 35px #233c3910}
*{box-sizing:border-box}body{margin:0;background:radial-gradient(ellipse at 82% 4%,#e0efe8 0,transparent 29rem),var(--paper);color:var(--ink);font:16px/1.5 Georgia,"Times New Roman",serif}
.shell{width:min(1120px,calc(100% - 40px));margin:auto}.mast{display:flex;justify-content:space-between;align-items:center;padding:22px 0;border-bottom:1px solid var(--line);font:700 12px/1.2 ui-monospace,monospace;letter-spacing:0;color:var(--muted);text-transform:uppercase}.mast strong{color:var(--green)}
header{padding:42px 0 30px}.eyebrow{font:700 12px/1.2 ui-monospace,monospace;letter-spacing:0;color:var(--green);text-transform:uppercase}h1{font-size:clamp(34px,5vw,58px);line-height:1.04;max-width:760px;margin:12px 0 14px;letter-spacing:0}header p{max-width:720px;color:var(--muted);font-size:19px;margin:0}.verdict{margin:8px 0 28px;padding:17px 20px;border-left:5px solid var(--lime);background:#fff;box-shadow:var(--shadow);font-size:18px}.verdict b{color:var(--green)}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:14px}.panel{background:var(--white);border:1px solid var(--line);border-radius:6px;padding:20px;box-shadow:var(--shadow)}.wide{grid-column:1/-1}.panel h2{font-size:23px;margin:0 0 8px;line-height:1.2}.panel p{margin:6px 0;color:var(--muted)}.toolbar{display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:10px;margin:14px 0 4px}.controls{display:flex;gap:5px}.controls button{font:700 12px ui-monospace,monospace;border:1px solid var(--line);background:#fff;color:var(--ink);padding:8px 11px;cursor:pointer;border-radius:4px}.controls button[aria-pressed=true]{background:var(--ink);border-color:var(--ink);color:#fff}.chart{width:100%;height:auto;overflow:visible}.chart text{font:12px ui-monospace,monospace;fill:var(--muted)}.chart .lang{font-weight:700;fill:var(--ink)}.chart .axis{stroke:#aab8b1;stroke-width:1}.chart .gridline{stroke:#e7ece8;stroke-width:1}.chart .range{stroke-width:2}.chart .box{stroke-width:1.5}.chart .median{stroke:var(--ink);stroke-width:2}.chart .sample{stroke:#fff;stroke-width:1}.key{display:flex;gap:18px;flex-wrap:wrap;font:12px ui-monospace,monospace;color:var(--muted);padding:0 0 8px}.key i{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:5px}.tablewrap{overflow-x:auto}table{width:100%;border-collapse:collapse;font:13px/1.35 ui-monospace,monospace;white-space:nowrap}th,td{text-align:right;padding:10px 8px;border-bottom:1px solid #e8ede8}th{color:var(--muted);font-size:11px;text-transform:uppercase}th:first-child,td:first-child{text-align:left}td:first-child{font-weight:700}.language{display:inline-flex;align-items:center;gap:8px}.swatch{width:10px;height:10px;border-radius:50%;display:inline-block}
.notes{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin:18px 0}.note{padding:16px 17px;border-radius:5px;background:#e8f1ec}.note:nth-child(2){background:#fff0e8}.note:nth-child(3){background:#eaf0fb}.note h3{margin:0 0 5px;font-size:17px}.note p{margin:0;color:#435654}.short{font:700 13px/1.45 ui-monospace,monospace;color:var(--ink)!important}.details{margin-top:12px}.details summary{font:700 13px ui-monospace,monospace;cursor:pointer;color:var(--green)}.raw{font:12px/1.8 ui-monospace,monospace;color:var(--muted);overflow-wrap:anywhere}.meta{display:grid;grid-template-columns:repeat(2,1fr);gap:9px 24px;margin-top:12px}.meta div{border-top:1px solid var(--line);padding-top:8px}.meta dt{font:700 11px ui-monospace,monospace;color:var(--muted);text-transform:uppercase}.meta dd{margin:3px 0 0;font:13px/1.4 ui-monospace,monospace;overflow-wrap:anywhere}.foot{padding:26px 0 45px;color:var(--muted);font:12px/1.6 ui-monospace,monospace}.small{font-size:13px!important}
@media(max-width:720px){.shell{width:min(100% - 24px,1120px)}.mast{padding:16px 0}header{padding:30px 0 20px}.grid{grid-template-columns:1fr}.wide{grid-column:auto}.notes{grid-template-columns:1fr}.panel{padding:15px}.meta{grid-template-columns:1fr}.toolbar{align-items:flex-start}.chart{min-width:560px}.chartframe{overflow-x:auto}.verdict{font-size:16px}}
.mast{border:0;color:#6b6d59}.mast strong{color:#dc4f37}.mast span:last-child{font-variant-numeric:tabular-nums}
.hero{position:relative;display:grid;grid-template-columns:minmax(0,1.5fr) minmax(220px,.7fr);gap:24px;align-items:center;margin:20px 0 16px;padding:32px 34px 34px;overflow:hidden;background:#d9ee77;border:1px solid #afc83f;border-radius:8px;box-shadow:0 12px 30px #4f5c2518;animation:flag-drop .55s cubic-bezier(.2,.8,.2,1) both}
.hero:after{content:"";position:absolute;left:0;right:0;bottom:0;height:12px;background:repeating-linear-gradient(135deg,#172a2b 0 12px,#f7f8f4 12px 24px)}
.hero-copy{position:relative;z-index:1}.hero .eyebrow{display:inline-flex;align-items:center;gap:8px;color:#294c43}.hero .eyebrow:before{content:"";width:9px;height:9px;border-radius:50%;background:#df4f37;box-shadow:0 0 0 4px #df4f3730}
.hero h1{max-width:640px;margin:12px 0 12px;font:800 clamp(42px,6.6vw,76px)/.93 "Trebuchet MS",Georgia,serif;letter-spacing:0;color:#172a2b}.hero p{max-width:620px;margin:0;color:#304c43;font-size:18px}.hero-tags{display:flex;flex-wrap:wrap;gap:8px;margin-top:20px}.hero-tags span{padding:6px 10px;border:1px solid #758f33;border-radius:99px;background:#eff8be;color:#304c43;font:700 11px/1.2 ui-monospace,monospace}
.hero-result{position:relative;z-index:1;justify-self:end;width:min(100%,280px);padding:17px 19px 18px;background:#172a2b;color:#f7f8f4;border-radius:6px;transform:rotate(1.2deg);box-shadow:7px 7px 0 #f7f8f4aa;animation:clock-pop .7s .2s cubic-bezier(.2,.8,.2,1) both}.hero-result .kicker{display:block;color:#d9ee77;font:700 10px ui-monospace,monospace;text-transform:uppercase}.hero-result strong{display:block;margin:5px 0;font:800 42px/.98 "Trebuchet MS",Georgia,sans-serif;font-variant-numeric:tabular-nums;letter-spacing:0}.hero-result .runner{font:700 13px ui-monospace,monospace}.hero-result .delta{display:block;margin-top:9px;color:#c5d3c8;font:11px/1.4 ui-monospace,monospace}
.race-call{margin:0 0 22px;padding:12px 15px;border-left:4px solid #df4f37;background:#fff0e8;color:#573a2b;font:700 13px/1.5 ui-monospace,monospace}.race-call b{color:#bd3b28}
.panel{box-shadow:0 9px 24px #263e3510;transition:transform .18s ease,box-shadow .18s ease}.panel:hover{transform:translateY(-2px);box-shadow:0 14px 28px #263e3518}.panel h2{font:800 24px/1.15 "Trebuchet MS",Georgia,serif}.panel h2:before{content:"// ";color:#df4f37;font:700 17px ui-monospace,monospace}
.race-board{display:grid;grid-template-columns:1fr;gap:7px;margin:14px 0 4px}.lane{display:grid;grid-template-columns:minmax(112px,160px) minmax(80px,1fr) 92px;gap:12px;align-items:center;padding:9px 11px;background:#f2f5ee;border:1px solid #e0e7dc;border-radius:4px;animation:lane-slide .42s both;animation-delay:calc(var(--lane) * 90ms)}.lane:first-child{background:#fff1cf;border-color:#e9c868}.lane-name{display:flex;align-items:center;gap:8px;font:800 13px ui-monospace,monospace}.lane-rank{display:grid;place-items:center;flex:0 0 23px;width:23px;height:23px;border-radius:50%;background:#172a2b;color:white;font-size:11px}.lane:first-child .lane-rank{background:#df4f37}.lane-track{height:11px;overflow:hidden;background:#dde5de;border-radius:99px}.lane-fill{height:100%;width:var(--bar);background:var(--lane-color);border-radius:99px;transform-origin:left;animation:bar-grow .8s cubic-bezier(.2,.8,.2,1) both;animation-delay:calc(var(--lane) * 100ms + 140ms)}.lane-time{text-align:right;font:800 13px ui-monospace,monospace;font-variant-numeric:tabular-nums}.lane-delta{grid-column:2/4;margin-top:-8px;color:#64706a;font:10px ui-monospace,monospace}
.grid{grid-template-columns:minmax(0,1fr) minmax(0,1fr)}.panel{min-width:0}.lane{min-width:0}.lane-track{min-width:0}
.key i{width:11px;height:11px}.chart .sample{stroke-width:1.5}.chart .gridline{stroke:#dbe4dc}.controls button{transition:background .15s ease,color .15s ease,transform .15s ease}.controls button:hover{transform:translateY(-1px)}
@keyframes flag-drop{from{opacity:0;transform:translateY(-12px)}to{opacity:1;transform:translateY(0)}}@keyframes clock-pop{0%{opacity:0;transform:translate(8px,6px) rotate(5deg) scale(.93)}100%{opacity:1;transform:rotate(1.2deg) scale(1)}}@keyframes lane-slide{from{opacity:0;transform:translateX(-10px)}to{opacity:1;transform:translateX(0)}}@keyframes bar-grow{from{transform:scaleX(0)}to{transform:scaleX(1)}}
@media(max-width:720px){.grid{grid-template-columns:minmax(0,1fr)}.hero{grid-template-columns:1fr;padding:24px 20px 30px;gap:18px;margin-top:10px}.hero h1{font-size:48px}.hero-result{justify-self:start;width:min(100%,250px)}.race-call{margin-bottom:16px}.lane{grid-template-columns:minmax(105px,125px) minmax(0,1fr) 76px;gap:8px;padding:8px}.lane-name,.lane-time{font-size:11px}.lane-delta{font-size:9px}}
@media(prefers-reduced-motion:reduce){*,*:before,*:after{animation-duration:.01ms!important;animation-iteration-count:1!important;scroll-behavior:auto!important;transition:none!important}}
</style>
</head>
<body>
<div class="shell">
<div class="mast"><span><strong>Synaptus rs2</strong> / baseline PSM benchmark</span><span id="runDate"></span></div>
<header class="hero"><div class="hero-copy"><div class="eyebrow">WireTargets · discrete PSM · local benchmark, not a world record</div><h1>Echoes in.<br>Stopwatch on.</h1><p>Four implementations process the same ultrasound scan. Both Rust ports are timed with the selected Rayon worker count; this covers baseline PSM, not NUFFT-PSM or SQ-PSM.</p><div class="hero-tags"><span id="lapTag"></span><span id="pixelTag"></span><span>same input · same math</span></div></div><aside class="hero-result" id="heroResult" aria-live="polite"></aside></header>
<div class="race-call" id="verdict" role="status"></div>
<main class="grid">
<section class="panel wide"><h2>Median board</h2><p>One-machine standings, sorted by median time. Bar length tracks elapsed time against the slowest median in this run; it is not a language score.</p><div class="race-board" id="raceBoard"></div></section>
<section class="panel wide"><h2>How long did each take?</h2><p>Every dot is a timed run. The thick bar is the middle half; the dark notch is the median. Lower means quicker.</p>
<div class="toolbar"><div class="controls" aria-label="Chart measure"><button data-mode="seconds" aria-pressed="true">Seconds</button><button data-mode="ratio" aria-pressed="false">Compared with rs2</button></div><div class="controls" aria-label="Axis scale"><button data-scale="linear" aria-pressed="true">Linear</button><button data-scale="log" aria-pressed="false">Log</button></div></div>
<div class="key" id="key"></div><div class="chartframe"><svg id="chart" class="chart" role="img" aria-label="Repeated benchmark run distribution"></svg></div></section>
<section class="panel wide"><h2>The numbers</h2><p>Q1 and Q3 bracket the middle 50% of runs. SD is the sample standard deviation. Quartiles use linear interpolation.</p><div class="tablewrap"><table><thead><tr><th>Implementation</th><th>Mean</th><th>Median</th><th>Q1</th><th>Q3</th><th>SD</th><th>Min</th><th>Max</th></tr></thead><tbody id="stats"></tbody></table></div>
<details class="details"><summary>Show every timed run</summary><div id="raw" class="raw"></div></details></section>
<section class="panel wide"><h2>Conclusions</h2><p id="conclusion"></p></section>
<section class="panel"><h2>What went on the clock?</h2><p>One end-to-end 2D PSM run on the bundled WireTargets dataset. Loading files and starting the programs are not timed.</p><dl class="meta" id="workload"></dl></section>
<section class="panel"><h2>Can we trust the race?</h2><p>All four complex output images are compared before timing results are reported. The error below is the largest pixel difference divided by rs2 Rust's peak pixel amplitude.</p><dl class="meta" id="agreement"></dl></section>
<section class="panel wide"><h2>Machine and tools</h2><dl class="meta" id="environment"></dl></section>
</main>
<footer class="foot">Benchmark tooling and report by GitHub Copilot, guided by Tomas. The latest result is benchmarks/benchmark_results.json; each run is archived at <code id="runRecord"></code>. Input and runner SHA-256 provenance is included. Rebuild with <code id="rebuild"></code>.</footer>
</div>
<script>
const REPORT = __REPORT_DATA__;
const keys = Object.keys(REPORT.results);
const palette = REPORT.palette;
const engineKey = REPORT.engine_key || (REPORT.environment.octave ? 'octave' : 'matlab');
const engineName = REPORT.engine_name || REPORT.languages[engineKey].name;
const fmt = value => `${value.toFixed(3)} s`;
document.getElementById('runDate').textContent = new Date(REPORT.timestamp).toLocaleString();
const rust = REPORT.results.rust.median_sec;
const fastest = keys.reduce((best, key) => REPORT.results[key].median_sec < REPORT.results[best].median_sec ? key : best, keys[0]);
const fastestName = REPORT.languages[fastest].name;
const fastestMedian = REPORT.results[fastest].median_sec;
const standings = [...keys].sort((a, b) => REPORT.results[a].median_sec - REPORT.results[b].median_sec);
const nextMedian = REPORT.results[standings[1]].median_sec;
const nextGap = nextMedian - fastestMedian;
const nextGapPercent = (nextGap / fastestMedian) * 100;
document.getElementById('lapTag').textContent = `${REPORT.warmups} warmups · ${REPORT.repeats} timed laps each`;
document.getElementById('pixelTag').textContent = `${REPORT.output_points.toLocaleString()} image pixels checked`;
document.getElementById('heroResult').innerHTML = `<span class="kicker">Quickest median · this machine</span><strong>${fastestMedian.toFixed(3)}<small> s</small></strong><span class="runner">${fastestName} · ${REPORT.repeats} timed laps</span><span class="delta">Next median: +${nextGap.toFixed(3)} s (${nextGapPercent.toFixed(1)}% longer)</span>`;
const worstError = Math.max(...Object.values(REPORT.agreement).map(value => value.max_abs_error_over_rs2_peak));
document.getElementById('verdict').innerHTML = `<b>Same-picture check passed.</b> Worst pixel difference: ${worstError.toExponential(2)} of rs2 Rust's peak. ${fastestName} had the shortest median in this run; one scan on one machine does not crown a universal champion.`;
const rs2 = REPORT.results.rust;
const rs1 = REPORT.results.rust_v1;
const rustGap = Math.abs(rs1.median_sec - rs2.median_sec) / Math.min(rs1.median_sec, rs2.median_sec) * 100;
const rustRangesOverlap = Math.max(rs1.q1_sec, rs2.q1_sec) <= Math.min(rs1.q3_sec, rs2.q3_sec);
const rustTake = rustRangesOverlap
    ? `Their middle-half ranges overlap, so the ${rustGap.toFixed(1)}% gap is a photo finish, not a coronation.`
    : `The observed middle-half ranges do not overlap; this is still one workload on one machine.`;
document.getElementById('conclusion').textContent = `With ${REPORT.rayon_threads} Rayon workers, rs2 clocked ${fmt(rs2.median_sec)} and rs clocked ${fmt(rs1.median_sec)}. ${rustTake} Octave and NumPy used one math-library thread, so this ranks these configurations, not languages. All four images agree within ${worstError.toExponential(2)} of rs2's peak across ${REPORT.output_points.toLocaleString()} pixels. The pixels have settled the argument; the stopwatch is still taking appeals.`;
const slowestMedian = Math.max(...keys.map(key => REPORT.results[key].median_sec));
document.getElementById('raceBoard').innerHTML = standings.map((key, index) => {
    const median = REPORT.results[key].median_sec;
    const gap = median - fastestMedian;
    const gapText = gap === 0 ? 'fastest median · local result' : `+${gap.toFixed(3)} s · +${(gap / fastestMedian * 100).toFixed(1)}% vs fastest`;
    return `<div class="lane" style="--lane:${index};--lane-color:${palette[key]}"><span class="lane-name"><i class="lane-rank">${index + 1}</i>${REPORT.languages[key].name}</span><div class="lane-track" role="img" aria-label="${REPORT.languages[key].name}: median ${fmt(median)}"><div class="lane-fill" style="--bar:${median / slowestMedian * 100}%"></div></div><span class="lane-time">${fmt(median)}</span><span class="lane-delta">${gapText}</span></div>`;
}).join('');
document.getElementById('key').innerHTML = keys.map(key => `<span><i style="background:${palette[key]}"></i>${REPORT.languages[key].name}</span>`).join('');
document.getElementById('stats').innerHTML = keys.map(key => {
  const s = REPORT.results[key];
  return `<tr><td><span class="language"><i class="swatch" style="background:${palette[key]}"></i>${REPORT.languages[key].name}</span></td><td>${fmt(s.mean_sec)}</td><td>${fmt(s.median_sec)}</td><td>${fmt(s.q1_sec)}</td><td>${fmt(s.q3_sec)}</td><td>${fmt(s.sd_sec)}</td><td>${fmt(s.min_sec)}</td><td>${fmt(s.max_sec)}</td></tr>`;
}).join('');
document.getElementById('raw').innerHTML = keys.map(key => `<p><b>${REPORT.languages[key].name}:</b> ${REPORT.results[key].samples_sec.map(fmt).join(' · ')}</p>`).join('');
const pairList = values => Object.entries(values).map(([label, value]) => `<div><dt>${label}</dt><dd>${value}</dd></div>`).join('');
document.getElementById('workload').innerHTML = pairList({
  'Input samples': `${REPORT.workload.nt.toLocaleString()} × ${REPORT.workload.nx.toLocaleString()}`,
  'Output pixels': REPORT.output_points.toLocaleString(),
  'Frequency band': `${REPORT.workload.f_low_hz / 1e6}–${REPORT.workload.f_high_hz / 1e6} MHz`,
  'Timed runs': `${REPORT.repeats} per implementation`,
  'Practice runs': `${REPORT.warmups} per implementation`,
  'Timer scope': 'FFT + migration + image reconstruction',
  'Excluded': 'File reads, process startup, report generation',
    'Rust builds': `Release, ${REPORT.rayon_threads} Rayon workers each`
});
document.getElementById('agreement').innerHTML = pairList(Object.fromEntries(keys.map(key => [
  REPORT.languages[key].name,
    `${REPORT.agreement[key].max_abs_error_over_rs2_peak.toExponential(2)} · ${REPORT.agreement[key].points.toLocaleString()} pixels`
])));
document.getElementById('environment').innerHTML = pairList({
  'CPU': REPORT.environment.cpu,
  'System': REPORT.environment.os,
    'Rust (both)': REPORT.environment.rustc,
    [engineName]: REPORT.environment[`${engineKey}_version`] || REPORT.environment[engineKey],
  'Python / NumPy': `${REPORT.environment.python} / ${REPORT.environment.numpy}`,
  'Parallelism': REPORT.environment.parallelism,
  'Git revision': REPORT.environment.git_revision
});
document.getElementById('rebuild').textContent = REPORT.rebuild_command || `python -I benchmarks/run_benchmark.py --engine ${engineKey}`;
document.getElementById('runRecord').textContent = REPORT.run_record;
let mode = 'seconds', scale = 'linear', dots = true;
function renderChart() {
  const svg = document.getElementById('chart');
  const width = 940, left = 148, right = 28, top = 40, row = 75;
  const height = top + row * keys.length + 36;
  const plotWidth = width - left - right;
  const values = keys.flatMap(key => REPORT.results[key].samples_sec.map(value => mode === 'seconds' ? value : value / rust));
  if (mode === 'ratio') values.push(1);
  let low = Math.min(...values), high = Math.max(...values);
  if (scale === 'linear') { const pad = (high - low || high * 0.1 || 1) * 0.12; low = Math.max(0, low - pad); high += pad; }
  else { low = Math.max(Number.MIN_VALUE, low * 0.78); high *= 1.28; }
  const project = value => {
    const fraction = scale === 'linear' ? (value - low) / (high - low) : (Math.log(value) - Math.log(low)) / (Math.log(high) - Math.log(low));
    return left + fraction * plotWidth;
  };
  const tickValues = Array.from({length:5}, (_, i) => scale === 'linear' ? low + (high - low) * i / 4 : Math.exp(Math.log(low) + (Math.log(high) - Math.log(low)) * i / 4));
  let content = '';
  tickValues.forEach((value, i) => {
    const x = project(value);
    content += `<line class="gridline" x1="${x}" y1="${top - 12}" x2="${x}" y2="${height - 28}"/><text x="${x}" y="${height - 7}" text-anchor="middle">${mode === 'seconds' ? `${value.toPrecision(2)} s` : `${value.toPrecision(2)}×`}</text>`;
  });
  if (mode === 'ratio' && 1 >= low && 1 <= high) content += `<line x1="${project(1)}" y1="${top - 12}" x2="${project(1)}" y2="${height - 28}" stroke="#172a2b" stroke-dasharray="4 4"/>`;
  keys.forEach((key, index) => {
    const s = REPORT.results[key];
    const vals = s.samples_sec.map(value => mode === 'seconds' ? value : value / rust);
    const y = top + index * row + 19;
    const q1 = mode === 'seconds' ? s.q1_sec : s.q1_sec / rust;
    const q3 = mode === 'seconds' ? s.q3_sec : s.q3_sec / rust;
    const median = mode === 'seconds' ? s.median_sec : s.median_sec / rust;
    const min = Math.min(...vals), max = Math.max(...vals), color = palette[key];
    content += `<text class="lang" x="0" y="${y + 4}">${REPORT.languages[key].name}</text>`;
    content += `<line class="range" x1="${project(min)}" y1="${y}" x2="${project(max)}" y2="${y}" stroke="${color}"/><line x1="${project(min)}" y1="${y-6}" x2="${project(min)}" y2="${y+6}" stroke="${color}"/><line x1="${project(max)}" y1="${y-6}" x2="${project(max)}" y2="${y+6}" stroke="${color}"/>`;
    content += `<rect class="box" x="${project(q1)}" y="${y-12}" width="${Math.max(1, project(q3)-project(q1))}" height="24" rx="3" fill="${color}" fill-opacity=".2" stroke="${color}"/><line class="median" x1="${project(median)}" y1="${y-12}" x2="${project(median)}" y2="${y+12}"/>`;
    if (dots) vals.forEach((value, n) => { const jitter = ((n % 5) - 2) * 3; content += `<circle class="sample" cx="${project(value)}" cy="${y+jitter}" r="3.5" fill="${color}"><title>Run ${n+1}: ${mode === 'seconds' ? fmt(value) : `${value.toFixed(3)}× rs2 median`}</title></circle>`; });
  });
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.innerHTML = content;
}
document.querySelectorAll('[data-mode]').forEach(button => button.addEventListener('click', () => {
  mode = button.dataset.mode;
  document.querySelectorAll('[data-mode]').forEach(item => item.setAttribute('aria-pressed', item === button));
  renderChart();
}));
document.querySelectorAll('[data-scale]').forEach(button => button.addEventListener('click', () => {
  scale = button.dataset.scale;
  document.querySelectorAll('[data-scale]').forEach(item => item.setAttribute('aria-pressed', item === button));
  renderChart();
}));
renderChart();
</script>
</body>
</html>'''
    return html.replace("__REPORT_DATA__", data)


def main():
    parser = argparse.ArgumentParser(description="Benchmark baseline WireTargets PSM in rs, rs2, MATLAB or Octave, and Python/NumPy.")
    parser.add_argument("--engine", choices=("matlab", "octave"), default="octave", help="third runtime to benchmark (default: octave)")
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=9)
    parser.add_argument("--rayon-threads", type=int, default=os.cpu_count() or 1, help="Rayon workers per Rust implementation (default: logical CPU count)")
    parser.add_argument("--render-only", action="store_true", help="rebuild the latest HTML from benchmark_results.json without timing workloads")
    args = parser.parse_args()
    if args.render_only:
        report = json.loads((OUTPUT / "benchmark_results.json").read_text())
        (OUTPUT / "performance_report.html").write_text(build_report(report))
        print(f"Rendered {OUTPUT / 'performance_report.html'} from the existing benchmark results.")
        return
    if args.warmups < 0 or args.repeats < 4 or args.rayon_threads < 1:
        parser.error("warmups must be non-negative, repeats at least 4, and Rayon threads at least 1")

    print(f"Running {args.warmups} warmups and {args.repeats} timed runs per implementation; Rayon={args.rayon_threads} for both Rust ports.", flush=True)
    with tempfile.TemporaryDirectory(prefix="synaptus-speed-") as temporary:
        temp_dir = Path(temporary)
        samples = run_drivers(args.engine, args.warmups, args.repeats, temp_dir, args.rayon_threads)
        agreement, output_points = check_images(args.engine, temp_dir)
    print("Image agreement passed; calculating distributions and writing the report.", flush=True)

    params = {}
    for line in (ROOT / "data/wire.txt").read_text().splitlines():
        parts = line.split()
        if parts:
            params[parts[0]] = [float(value) for value in parts[1:]]
    record_name = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S-%f") + "-psm-benchmark"
    record_path = Path("runs") / record_name
    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "run_record": record_path.as_posix(),
        "benchmark_scope": "Discrete baseline PSM parity across synaptus-rs, synaptus-rs2, and the selected third runtime; NUFFT-PSM and SQ-PSM are not timed.",
        "engine_key": args.engine,
        "engine_name": LANGUAGES[args.engine][0],
        "rebuild_command": f".venv/bin/python -I benchmarks/run_benchmark.py --engine {args.engine} --warmups {args.warmups} --repeats {args.repeats} --rayon-threads {args.rayon_threads}",
        "warmups": args.warmups,
        "repeats": args.repeats,
        "rayon_threads": args.rayon_threads,
        "tolerance_normalized_to_rs2_peak": TOLERANCE,
        "palette": {key: LANGUAGES[key][1] for key in samples},
        "languages": {key: {"name": LANGUAGES[key][0]} for key in samples},
        "results": {key: describe(values) for key, values in samples.items()},
        "agreement": agreement,
        "output_points": output_points,
        "workload": {
            "dataset": "WireTargets (bundled data/wire.f64)",
            "nt": int(params["dims"][0]),
            "nx": int(params["dims"][1]),
            "fs_hz": params["fs"][0],
            "f_low_hz": params["fLow"][0],
            "f_high_hz": params["fHigh"][0],
            "velocity_m_s": params["cc"][0],
            "thickness_m": params["thick"][0],
            "x_step_m": params["xStep"][0],
            "depth_step": "(c / 2) / f_high, matching PsmOptions default",
        },
        "environment": collect_environment(args.engine, args.rayon_threads),
        "provenance": {
            name: sha256(ROOT / path)
            for name, path in {
                "wire_f64_sha256": "data/wire.f64",
                "wire_metadata_sha256": "data/wire.txt",
                "rust_runner_sha256": "examples/speed_benchmark.rs",
                "rs1_rust_runner_sha256": "../synaptus-rs/examples/speed_benchmark.rs",
                "rs1_psm_implementation_sha256": "../synaptus-rs/src/psm.rs",
                "rs1_cargo_manifest_sha256": "../synaptus-rs/Cargo.toml",
                "rs1_cargo_lock_sha256": "../synaptus-rs/Cargo.lock",
                f"{args.engine}_runner_sha256": f"benchmarks/psm_benchmark{'_matlab' if args.engine == 'matlab' else ''}.m",
                "python_runner_sha256": "benchmarks/python_psm_benchmark.py",
                "orchestrator_sha256": "benchmarks/run_benchmark.py",
                "psm_implementation_sha256": "src/psm.rs",
                "cargo_manifest_sha256": "Cargo.toml",
                "cargo_lock_sha256": "Cargo.lock",
                "benchmark_requirements_sha256": "requirements-benchmark.txt",
            }.items()
        },
    }
    result_json = json.dumps(report, indent=2, allow_nan=False) + "\n"
    report_html = build_report(report)
    (OUTPUT / "benchmark_results.json").write_text(result_json)
    (OUTPUT / "performance_report.html").write_text(report_html)
    archived = ROOT / record_path
    archived.mkdir(parents=True, exist_ok=False)
    (archived / "benchmark_results.json").write_text(result_json)
    (archived / "performance_report.html").write_text(report_html)
    print(f"Wrote {OUTPUT / 'performance_report.html'}")
    print(f"Wrote {OUTPUT / 'benchmark_results.json'}")
    for key, values in samples.items():
        print(f"{LANGUAGES[key][0]} median: {statistics.median(values):.4f} s")


if __name__ == "__main__":
    main()