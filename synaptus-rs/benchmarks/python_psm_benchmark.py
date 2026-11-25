import argparse
import time
from pathlib import Path

import numpy as np


def read_params(path):
    params = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if parts:
            params[parts[0]] = np.array([float(value) for value in parts[1:]])
    return params


def psm(data, params):
    nt, nx = map(int, params["dims"][:2])
    fs, delay, velocity = params["fs"][0], params["tDelay"][0], params["cc"][0]
    f_low, f_high = params["fLow"][0], params["fHigh"][0]
    x_step, thickness = params["xStep"][0], params["thick"][0]
    nfft_t = 1 << (nt - 1).bit_length()
    nfft_x = 1 << (nx - 1).bit_length()

    time_index = np.arange(nfft_t)
    time_index[time_index >= (nfft_t + 1) // 2] -= nfft_t
    omega = time_index * (2.0 * np.pi * fs / nfft_t)
    keep = (omega >= 2.0 * np.pi * f_low) & (omega <= 2.0 * np.pi * f_high)
    omega = omega[keep]
    spectrum = np.fft.fft(data, n=nfft_t, axis=0)[keep, :]
    spectrum *= np.exp(-1j * omega * delay)[:, None]
    spectrum = np.fft.fft(spectrum, n=nfft_x, axis=1).T

    space_index = np.arange(nfft_x)
    space_index[space_index >= (nfft_x + 1) // 2] -= nfft_x
    kx = space_index * ((2.0 * np.pi / x_step) / nfft_x)
    kz2 = (2.0 / velocity) ** 2 * omega[None, :] ** 2 - kx[:, None] ** 2
    propagating = kz2 >= 0.0
    kz = np.sqrt(np.maximum(kz2, 0.0))
    spectrum[~propagating] = 0.0

    dz = (velocity / 2.0) / f_high
    z_offset = delay * (velocity / 2.0)
    nz = int(np.ceil((thickness - z_offset) / dz)) + 1
    wave = spectrum * np.exp(1j * kz * z_offset)
    phase = np.exp(1j * kz * dz)
    image = np.empty((nz, nx), dtype=np.complex128)
    for iz in range(nz):
        image[iz] = np.fft.ifft(np.sum(wave, axis=1))[:nx]
        wave *= phase
    return image


def save_image(path, image):
    values = image.reshape(-1)
    raw = np.empty(values.size * 2, dtype="<f8")
    raw[0::2], raw[1::2] = values.real, values.imag
    raw.tofile(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("warmups", type=int)
    parser.add_argument("repeats", type=int)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    params = read_params(args.data_dir / "wire.txt")
    nt, nx = map(int, params["dims"][:2])
    raw = np.fromfile(args.data_dir / "wire.f64", dtype="<f8")
    data = raw.reshape((nt, nx), order="F")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for _ in range(args.warmups):
        psm(data, params)
    samples = []
    result = None
    for _ in range(args.repeats):
        started = time.perf_counter()
        result = psm(data, params)
        samples.append(time.perf_counter() - started)
    save_image(args.output_dir / "python_result.c64", result)
    (args.output_dir / "python_samples.csv").write_text(",".join(f"{value:.9f}" for value in samples) + "\n")
    print(f"Python/NumPy PSM: {len(samples)} samples")


if __name__ == "__main__":
    main()