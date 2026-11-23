% Check: MULOK's band upsampling Poxy = fft(ifft(Poxy), up*n) and the omega
% grid assigned to it in mulok.m / learn_stolt.m.
% Oracle: the band samples a known analytic spectrum X(w) = exp(-i*w*tau)
% (a pure delay tau). Interpolated values must equal X at the TRUE grid.
fs = 12.5e6; nFFTt = 1024; up = 4; tau = 3e-6;
dOmega = 2*pi*fs/nFFTt;
omega = ifftshift(((0:nFFTt-1) - floor(nFFTt/2))' * dOmega);
band = (omega >= -2*pi*4e6) & (omega <= -2*pi*1e6);   % mulok.m convention
omegaBand = omega(band); n = numel(omegaBand);
X = exp(-1i*omegaBand*tau);                           % analytic band samples
Xup = fft(ifft(X), up*n);                             % mulok.m line
tmp = (0:(n*up-1))*(dOmega/up);
grid_code = -tmp(end:-1:1).' + omegaBand(end);        % mulok.m grid
grid_true = omegaBand(1) + (0:(n*up-1)).'*(dOmega/up);% DFT-interpolation grid
printf('grid offset code-true (in units of dOmega): %.6f (expected -(1-1/up) = %.6f)\n', ...
  mean(grid_code-grid_true)/dOmega, -(1-1/up));
inner = 1:(up*(n-1)+1);   % samples not affected by wrap-around between band ends
e_true = max(abs(Xup(inner) - exp(-1i*grid_true(inner)*tau)));
e_code = max(abs(Xup(inner) - exp(-1i*grid_code(inner)*tau)));
printf('max |Xup - X(grid_true)| = %.3e   max |Xup - X(grid_code)| = %.3e\n', e_true, e_code);
printf('exact-sample check (every up-th sample equals input): %.3e\n', max(abs(Xup(1:up:end)-X)));
