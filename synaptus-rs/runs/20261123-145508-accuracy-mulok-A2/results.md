dataset | configuration | complex raw | complex fit | envelope raw | envelope fit
---|---|---|---|---|---
synth1 | PSM self-check: nfft_t x2 (scaled 1/2) vs x1 | 3.656e-6 | 3.656e-6 | 3.424e-6 | 3.424e-6
synth1 | legacy mulok.m (vs PSM, its dz grid) | 1.001e0 | 9.131e-1 | 9.949e-1 | 1.411e-2
synth1 | accurate baseline (cubic, up=4, demod) | 3.696e-6 | 3.696e-6 | 3.582e-6 | 3.582e-6
synth1 | baseline + legacy squared Jacobian | 9.003e-3 | 7.408e-3 | 8.891e-3 | 7.270e-3
synth1 | baseline + legacy kz truncation | 3.696e-6 | 3.696e-6 | 3.583e-6 | 3.583e-6
synth1 | baseline + legacy -0.75 bin ω offset | 1.166e0 | 5.442e-2 | 5.484e-3 | 5.484e-3
synth1 | baseline + all three legacy defects | 1.164e0 | 5.530e-2 | 9.406e-3 | 7.967e-3
synth1 | thesis as coded in mulok.m interp (linear, up=4, no demod) | 2.409e-1 | 1.400e-1 | 2.409e-1 | 1.400e-1
synth2 | PSM self-check: nfft_t x2 (scaled 1/2) vs x1 | 3.623e-4 | 3.623e-4 | 3.093e-4 | 3.093e-4
synth2 | legacy mulok.m (vs PSM, its dz grid) | 9.971e-1 | 8.713e-1 | 9.937e-1 | 5.718e-2
synth2 | accurate baseline (cubic, up=4, demod) | 3.942e-4 | 3.928e-4 | 3.560e-4 | 3.545e-4
synth2 | baseline + legacy squared Jacobian | 3.221e-2 | 2.704e-2 | 2.934e-2 | 2.348e-2
synth2 | baseline + legacy kz truncation | 4.849e-4 | 4.838e-4 | 4.023e-4 | 4.009e-4
synth2 | baseline + legacy -0.75 bin ω offset | 1.004e0 | 5.007e-2 | 1.608e-2 | 1.608e-2
synth2 | baseline + all three legacy defects | 9.952e-1 | 5.819e-2 | 3.202e-2 | 2.707e-2
synth2 | thesis as coded in mulok.m interp (linear, up=4, no demod) | 3.750e-2 | 3.643e-2 | 3.739e-2 | 3.632e-2
wire | PSM self-check: nfft_t x2 (scaled 1/2) vs x1 | 1.728e-2 | 1.728e-2 | 1.208e-2 | 1.208e-2
wire | legacy mulok.m (vs PSM, its dz grid) | 9.997e-1 | 9.699e-1 | 9.986e-1 | 4.947e-2
wire | accurate baseline (cubic, up=4, demod) | 2.198e-2 | 2.197e-2 | 1.526e-2 | 1.525e-2
wire | baseline + legacy squared Jacobian | 2.462e-2 | 2.429e-2 | 1.755e-2 | 1.705e-2
wire | baseline + legacy kz truncation | 2.392e-2 | 2.392e-2 | 1.661e-2 | 1.660e-2
wire | baseline + legacy -0.75 bin ω offset | 1.478e0 | 8.636e-1 | 2.243e-2 | 2.243e-2
wire | baseline + all three legacy defects | 1.475e0 | 8.636e-1 | 2.280e-2 | 2.242e-2
wire | thesis as coded in mulok.m interp (linear, up=4, no demod) | 6.920e-2 | 5.599e-2 | 6.495e-2 | 5.041e-2
pins | PSM self-check: nfft_t x2 (scaled 1/2) vs x1 | 2.172e-3 | 2.172e-3 | 1.527e-3 | 1.527e-3
pins | legacy mulok.m (vs PSM, its dz grid) | 9.999e-1 | 9.989e-1 | 9.980e-1 | 3.046e-2
pins | accurate baseline (cubic, up=4, demod) | 4.136e-3 | 4.135e-3 | 2.750e-3 | 2.747e-3
pins | baseline + legacy squared Jacobian | 8.015e-3 | 7.951e-3 | 5.979e-3 | 5.890e-3
pins | baseline + legacy kz truncation | 1.108e-2 | 1.108e-2 | 7.723e-3 | 7.721e-3
pins | baseline + legacy -0.75 bin ω offset | 7.561e-1 | 6.865e-1 | 6.243e-3 | 6.241e-3
pins | baseline + all three legacy defects | 7.561e-1 | 6.868e-1 | 1.066e-2 | 1.061e-2
pins | thesis as coded in mulok.m interp (linear, up=4, no demod) | 3.336e-2 | 3.208e-2 | 3.142e-2 | 3.004e-2
plane | PSM self-check: nfft_t x2 (scaled 1/2) vs x1 | 6.228e-3 | 6.228e-3 | 4.197e-3 | 4.197e-3
plane | legacy mulok.m (vs PSM, its dz grid) | 9.995e-1 | 9.887e-1 | 9.960e-1 | 5.374e-2
plane | accurate baseline (cubic, up=4, demod) | 6.441e-3 | 6.435e-3 | 4.453e-3 | 4.446e-3
plane | baseline + legacy squared Jacobian | 1.215e-2 | 1.202e-2 | 9.851e-3 | 9.698e-3
plane | baseline + legacy kz truncation | 1.659e-2 | 1.659e-2 | 1.158e-2 | 1.158e-2
plane | baseline + legacy -0.75 bin ω offset | 4.665e-1 | 4.324e-1 | 1.029e-2 | 1.028e-2
plane | baseline + all three legacy defects | 4.659e-1 | 4.322e-1 | 1.630e-2 | 1.619e-2
plane | thesis as coded in mulok.m interp (linear, up=4, no demod) | 4.319e-2 | 4.238e-2 | 4.209e-2 | 4.125e-2
