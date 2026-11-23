dataset | configuration | raw rel.err | scale-fitted rel.err | time [s]
---|---|---|---|---
synth1 | legacy mulok.m (vs PSM on its dz grid) | 1.001e0 | 9.131e-1 | 0.01
synth1 | thesis: linear, up=4 | 2.409e-1 | 1.400e-1 | 0.01
synth1 | thesis + legacy squared Jacobian | 2.444e-1 | 1.403e-1 | 0.01
synth1 | thesis + legacy kz truncation | 2.409e-1 | 1.400e-1 | 0.01
synth1 | thesis + legacy -0.75 bin ω offset | 1.062e0 | 1.499e-1 | 0.01
synth1 | thesis + all three legacy defects | 1.060e0 | 1.504e-1 | 0.01
synth1 | interp=Linear up=1 demod=false | 1.400e0 | 1.000e0 | 0.00
synth1 | interp=Linear up=2 demod=false | 6.518e-1 | 1.784e-1 | 0.00
synth1 | interp=Linear up=4 demod=false | 2.409e-1 | 1.400e-1 | 0.01
synth1 | interp=Linear up=8 demod=false | 5.780e-2 | 7.621e-3 | 0.01
synth1 | interp=Cubic up=1 demod=false | 1.414e0 | 1.000e0 | 0.00
synth1 | interp=Cubic up=2 demod=false | 6.480e-1 | 1.495e-1 | 0.00
synth1 | interp=Cubic up=4 demod=false | 1.055e-1 | 5.716e-2 | 0.01
synth1 | interp=Cubic up=8 demod=false | 7.135e-3 | 8.841e-4 | 0.01
synth1 | interp=Lanczos(6) up=1 demod=false | 1.414e0 | 1.000e0 | 0.01
synth1 | interp=Lanczos(6) up=2 demod=false | 7.352e-1 | 2.267e-1 | 0.01
synth1 | interp=Lanczos(6) up=4 demod=false | 5.050e-3 | 3.637e-3 | 0.01
synth1 | interp=Lanczos(6) up=8 demod=false | 8.207e-4 | 1.151e-4 | 0.01
synth1 | interp=Linear up=1 demod=true | 1.994e-3 | 1.464e-3 | 0.00
synth1 | interp=Linear up=2 demod=true | 4.634e-4 | 3.167e-4 | 0.00
synth1 | interp=Linear up=4 demod=true | 1.245e-4 | 9.117e-5 | 0.01
synth1 | interp=Linear up=8 demod=true | 2.923e-5 | 2.012e-5 | 0.01
synth1 | interp=Cubic up=1 demod=true | 3.030e-5 | 2.939e-5 | 0.00
synth1 | interp=Cubic up=2 demod=true | 4.229e-6 | 4.204e-6 | 0.00
synth1 | interp=Cubic up=4 demod=true | 3.696e-6 | 3.696e-6 | 0.01
synth1 | interp=Cubic up=8 demod=true | 3.694e-6 | 3.694e-6 | 0.01
synth1 | interp=Lanczos(6) up=1 demod=true | 8.561e-4 | 7.733e-4 | 0.01
synth1 | interp=Lanczos(6) up=2 demod=true | 3.831e-4 | 7.582e-5 | 0.01
synth1 | interp=Lanczos(6) up=4 demod=true | 5.032e-4 | 3.342e-4 | 0.01
synth1 | interp=Lanczos(6) up=8 demod=true | 3.768e-4 | 2.278e-5 | 0.02
synth2 | legacy mulok.m (vs PSM on its dz grid) | 9.971e-1 | 8.713e-1 | 0.01
synth2 | thesis: linear, up=4 | 3.750e-2 | 3.643e-2 | 0.01
synth2 | thesis + legacy squared Jacobian | 5.293e-2 | 4.606e-2 | 0.01
synth2 | thesis + legacy kz truncation | 3.750e-2 | 3.644e-2 | 0.01
synth2 | thesis + legacy -0.75 bin ω offset | 1.000e0 | 6.187e-2 | 0.01
synth2 | thesis + all three legacy defects | 9.910e-1 | 6.906e-2 | 0.01
synth2 | interp=Linear up=1 demod=false | 1.537e-1 | 1.818e-2 | 0.00
synth2 | interp=Linear up=2 demod=false | 4.278e-1 | 4.056e-1 | 0.00
synth2 | interp=Linear up=4 demod=false | 3.750e-2 | 3.643e-2 | 0.01
synth2 | interp=Linear up=8 demod=false | 1.209e-2 | 1.185e-2 | 0.01
synth2 | interp=Cubic up=1 demod=false | 4.965e-2 | 9.503e-3 | 0.00
synth2 | interp=Cubic up=2 demod=false | 4.623e-1 | 4.209e-1 | 0.00
synth2 | interp=Cubic up=4 demod=false | 1.349e-2 | 1.349e-2 | 0.01
synth2 | interp=Cubic up=8 demod=false | 1.399e-3 | 1.399e-3 | 0.01
synth2 | interp=Lanczos(6) up=1 demod=false | 1.673e-3 | 1.186e-3 | 0.01
synth2 | interp=Lanczos(6) up=2 demod=false | 5.741e-1 | 4.983e-1 | 0.01
synth2 | interp=Lanczos(6) up=4 demod=false | 4.185e-4 | 3.978e-4 | 0.01
synth2 | interp=Lanczos(6) up=8 demod=false | 7.281e-4 | 6.344e-4 | 0.02
synth2 | interp=Linear up=1 demod=true | 6.223e-2 | 1.962e-2 | 0.00
synth2 | interp=Linear up=2 demod=true | 1.654e-2 | 6.767e-3 | 0.01
synth2 | interp=Linear up=4 demod=true | 4.076e-3 | 1.459e-3 | 0.01
synth2 | interp=Linear up=8 demod=true | 1.086e-3 | 5.160e-4 | 0.02
synth2 | interp=Cubic up=1 demod=true | 8.611e-3 | 3.671e-3 | 0.00
synth2 | interp=Cubic up=2 demod=true | 7.385e-4 | 5.277e-4 | 0.01
synth2 | interp=Cubic up=4 demod=true | 3.942e-4 | 3.928e-4 | 0.01
synth2 | interp=Cubic up=8 demod=true | 3.924e-4 | 3.924e-4 | 0.02
synth2 | interp=Lanczos(6) up=1 demod=true | 9.534e-4 | 6.763e-4 | 0.01
synth2 | interp=Lanczos(6) up=2 demod=true | 5.750e-4 | 5.674e-4 | 0.01
synth2 | interp=Lanczos(6) up=4 demod=true | 8.033e-4 | 7.304e-4 | 0.01
synth2 | interp=Lanczos(6) up=8 demod=true | 6.307e-4 | 5.082e-4 | 0.03
wire | legacy mulok.m (vs PSM on its dz grid) | 9.997e-1 | 9.699e-1 | 0.06
wire | thesis: linear, up=4 | 6.920e-2 | 5.599e-2 | 0.18
wire | thesis + legacy squared Jacobian | 7.219e-2 | 5.713e-2 | 0.23
wire | thesis + legacy kz truncation | 6.984e-2 | 5.678e-2 | 0.19
wire | thesis + legacy -0.75 bin ω offset | 1.430e0 | 8.543e-1 | 0.16
wire | thesis + all three legacy defects | 1.427e0 | 8.543e-1 | 0.20
wire | interp=Linear up=1 demod=false | 6.311e-1 | 6.296e-1 | 0.02
wire | interp=Linear up=2 demod=false | 2.197e-1 | 1.795e-1 | 0.03
wire | interp=Linear up=4 demod=false | 6.920e-2 | 5.599e-2 | 0.20
wire | interp=Linear up=8 demod=false | 2.722e-2 | 2.502e-2 | 0.58
wire | interp=Cubic up=1 demod=false | 6.614e-1 | 6.588e-1 | 0.02
wire | interp=Cubic up=2 demod=false | 1.491e-1 | 1.323e-1 | 0.04
wire | interp=Cubic up=4 demod=false | 2.669e-2 | 2.566e-2 | 0.23
wire | interp=Cubic up=8 demod=false | 2.199e-2 | 2.199e-2 | 0.58
wire | interp=Lanczos(6) up=1 demod=false | 8.047e-1 | 7.545e-1 | 0.05
wire | interp=Lanczos(6) up=2 demod=false | 4.042e-2 | 4.015e-2 | 0.06
wire | interp=Lanczos(6) up=4 demod=false | 2.203e-2 | 2.202e-2 | 0.22
wire | interp=Lanczos(6) up=8 demod=false | 2.200e-2 | 2.200e-2 | 0.60
wire | interp=Linear up=1 demod=true | 2.527e-1 | 1.877e-1 | 0.03
wire | interp=Linear up=2 demod=true | 7.847e-2 | 5.801e-2 | 0.04
wire | interp=Linear up=4 demod=true | 2.865e-2 | 2.504e-2 | 0.18
wire | interp=Linear up=8 demod=true | 2.233e-2 | 2.209e-2 | 0.57
wire | interp=Cubic up=1 demod=true | 1.583e-1 | 1.194e-1 | 0.02
wire | interp=Cubic up=2 demod=true | 2.628e-2 | 2.435e-2 | 0.04
wire | interp=Cubic up=4 demod=true | 2.198e-2 | 2.197e-2 | 0.21
wire | interp=Cubic up=8 demod=true | 2.199e-2 | 2.198e-2 | 0.58
wire | interp=Lanczos(6) up=1 demod=true | 2.310e-2 | 2.298e-2 | 0.05
wire | interp=Lanczos(6) up=2 demod=true | 2.207e-2 | 2.206e-2 | 0.06
wire | interp=Lanczos(6) up=4 demod=true | 2.199e-2 | 2.199e-2 | 0.21
wire | interp=Lanczos(6) up=8 demod=true | 2.200e-2 | 2.199e-2 | 0.61
pins | legacy mulok.m (vs PSM on its dz grid) | 9.999e-1 | 9.989e-1 | 0.02
pins | thesis: linear, up=4 | 3.336e-2 | 3.208e-2 | 0.01
pins | thesis + legacy squared Jacobian | 3.411e-2 | 3.259e-2 | 0.01
pins | thesis + legacy kz truncation | 3.491e-2 | 3.368e-2 | 0.01
pins | thesis + legacy -0.75 bin ω offset | 7.315e-1 | 6.701e-1 | 0.01
pins | thesis + all three legacy defects | 7.315e-1 | 6.704e-1 | 0.01
pins | interp=Linear up=1 demod=false | 6.520e-1 | 5.840e-1 | 0.00
pins | interp=Linear up=2 demod=false | 1.471e-1 | 1.463e-1 | 0.01
pins | interp=Linear up=4 demod=false | 3.336e-2 | 3.208e-2 | 0.01
pins | interp=Linear up=8 demod=false | 1.091e-2 | 1.061e-2 | 0.03
pins | interp=Cubic up=1 demod=false | 7.421e-1 | 6.341e-1 | 0.01
pins | interp=Cubic up=2 demod=false | 1.092e-1 | 1.090e-1 | 0.01
pins | interp=Cubic up=4 demod=false | 6.839e-3 | 6.668e-3 | 0.01
pins | interp=Cubic up=8 demod=false | 4.172e-3 | 4.169e-3 | 0.02
pins | interp=Lanczos(6) up=1 demod=false | 9.125e-1 | 7.191e-1 | 0.02
pins | interp=Lanczos(6) up=2 demod=false | 2.201e-2 | 2.201e-2 | 0.02
pins | interp=Lanczos(6) up=4 demod=false | 4.162e-3 | 4.156e-3 | 0.03
pins | interp=Lanczos(6) up=8 demod=false | 4.196e-3 | 4.189e-3 | 0.04
pins | interp=Linear up=1 demod=true | 6.757e-2 | 5.564e-2 | 0.01
pins | interp=Linear up=2 demod=true | 1.904e-2 | 1.588e-2 | 0.01
pins | interp=Linear up=4 demod=true | 6.010e-3 | 5.344e-3 | 0.01
pins | interp=Linear up=8 demod=true | 4.293e-3 | 4.229e-3 | 0.03
pins | interp=Cubic up=1 demod=true | 1.912e-2 | 1.700e-2 | 0.01
pins | interp=Cubic up=2 demod=true | 4.394e-3 | 4.339e-3 | 0.01
pins | interp=Cubic up=4 demod=true | 4.136e-3 | 4.135e-3 | 0.02
pins | interp=Cubic up=8 demod=true | 4.135e-3 | 4.134e-3 | 0.03
pins | interp=Lanczos(6) up=1 demod=true | 4.263e-3 | 4.263e-3 | 0.02
pins | interp=Lanczos(6) up=2 demod=true | 4.192e-3 | 4.192e-3 | 0.02
pins | interp=Lanczos(6) up=4 demod=true | 4.145e-3 | 4.137e-3 | 0.03
pins | interp=Lanczos(6) up=8 demod=true | 4.180e-3 | 4.168e-3 | 0.04
plane | legacy mulok.m (vs PSM on its dz grid) | 9.995e-1 | 9.887e-1 | 0.55
plane | thesis: linear, up=4 | 4.319e-2 | 4.238e-2 | 0.48
plane | thesis + legacy squared Jacobian | 4.472e-2 | 4.365e-2 | 0.43
plane | thesis + legacy kz truncation | 4.585e-2 | 4.509e-2 | 0.41
plane | thesis + legacy -0.75 bin ω offset | 4.519e-1 | 4.173e-1 | 0.40
plane | thesis + all three legacy defects | 4.514e-1 | 4.172e-1 | 0.44
plane | interp=Linear up=1 demod=false | 2.160e-1 | 2.151e-1 | 0.21
plane | interp=Linear up=2 demod=false | 1.438e-1 | 1.433e-1 | 0.28
plane | interp=Linear up=4 demod=false | 4.319e-2 | 4.238e-2 | 0.42
plane | interp=Linear up=8 demod=false | 1.326e-2 | 1.298e-2 | 0.75
plane | interp=Cubic up=1 demod=false | 2.095e-1 | 2.095e-1 | 0.22
plane | interp=Cubic up=2 demod=false | 1.357e-1 | 1.354e-1 | 0.28
plane | interp=Cubic up=4 demod=false | 1.894e-2 | 1.866e-2 | 0.40
plane | interp=Cubic up=8 demod=false | 6.609e-3 | 6.598e-3 | 0.65
plane | interp=Lanczos(6) up=1 demod=false | 2.081e-1 | 2.081e-1 | 0.60
plane | interp=Lanczos(6) up=2 demod=false | 1.419e-1 | 1.417e-1 | 0.67
plane | interp=Lanczos(6) up=4 demod=false | 6.515e-3 | 6.512e-3 | 0.79
plane | interp=Lanczos(6) up=8 demod=false | 6.478e-3 | 6.471e-3 | 1.02
plane | interp=Linear up=1 demod=true | 1.286e-1 | 7.669e-2 | 0.24
plane | interp=Linear up=2 demod=true | 3.718e-2 | 2.425e-2 | 0.30
plane | interp=Linear up=4 demod=true | 1.077e-2 | 7.903e-3 | 0.44
plane | interp=Linear up=8 demod=true | 6.840e-3 | 6.557e-3 | 0.76
plane | interp=Cubic up=1 demod=true | 4.233e-2 | 2.652e-2 | 0.24
plane | interp=Cubic up=2 demod=true | 7.309e-3 | 6.855e-3 | 0.31
plane | interp=Cubic up=4 demod=true | 6.441e-3 | 6.435e-3 | 0.46
plane | interp=Cubic up=8 demod=true | 6.435e-3 | 6.433e-3 | 0.78
plane | interp=Lanczos(6) up=1 demod=true | 6.546e-3 | 6.511e-3 | 0.70
plane | interp=Lanczos(6) up=2 demod=true | 6.591e-3 | 6.559e-3 | 0.77
plane | interp=Lanczos(6) up=4 demod=true | 6.439e-3 | 6.436e-3 | 0.89
plane | interp=Lanczos(6) up=8 demod=true | 6.474e-3 | 6.465e-3 | 1.19
