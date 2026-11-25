use std::path::Path;
use std::time::Instant;
use synaptus::psm::{psm, PsmOptions};
use synaptus::{io, FocusedImage};

fn save_image(path: &Path, image: &FocusedImage) {
    let mut bytes = Vec::with_capacity(image.layers[0].data.len() * 16);
    for value in &image.layers[0].data {
        bytes.extend_from_slice(&value.re.to_le_bytes());
        bytes.extend_from_slice(&value.im.to_le_bytes());
    }
    std::fs::write(path, bytes).expect("write result image");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 5, "usage: speed_benchmark DATA_DIR OUTPUT_DIR WARMUPS REPEATS");
    let data_dir = Path::new(&args[1]);
    let output_dir = Path::new(&args[2]);
    let warmups: usize = args[3].parse().expect("warmups must be an integer");
    let repeats: usize = args[4].parse().expect("repeats must be an integer");
    assert!(repeats >= 4, "use at least four timed repetitions");
    std::fs::create_dir_all(output_dir).expect("create output directory");
    let (data, acquisition, medium) = io::load_dataset(data_dir, "wire").expect("load WireTargets");

    for _ in 0..warmups {
        std::hint::black_box(psm(&data, &acquisition, &medium, &PsmOptions::default()).unwrap());
    }
    let mut samples = Vec::with_capacity(repeats);
    let mut last = None;
    for _ in 0..repeats {
        let started = Instant::now();
        let image = psm(&data, &acquisition, &medium, &PsmOptions::default()).unwrap();
        samples.push(started.elapsed().as_secs_f64());
        std::hint::black_box(image.layers[0].data[0].norm());
        last = Some(image);
    }
    save_image(&output_dir.join("rust_result.c64"), last.as_ref().unwrap());
    let samples = samples.iter().map(|x| format!("{x:.9}")).collect::<Vec<_>>().join(",");
    std::fs::write(output_dir.join("rust_samples.csv"), format!("{samples}\n"))
        .expect("write timing samples");
    println!("Rust PSM: {} samples", repeats);
}