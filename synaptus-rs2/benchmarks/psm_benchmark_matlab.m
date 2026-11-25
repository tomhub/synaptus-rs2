function psm_benchmark_matlab(data_dir, output_dir, warmups, repeats)
  if repeats < 1
    error('repeats must be positive');
  end
  if exist('maxNumCompThreads', 'builtin') || exist('maxNumCompThreads', 'file')
    try
      maxNumCompThreads(1);
    catch thread_error
      warning('Could not limit MATLAB worker threads: %s', thread_error.message);
    end
  end

  params = read_params(fullfile(data_dir, 'wire.txt'));
  nt = params.dims(1);
  nx = params.dims(2);
  fid = fopen(fullfile(data_dir, 'wire.f64'), 'rb', 'ieee-le');
  if fid < 0
    error('Cannot open WireTargets data.');
  end
  data = reshape(fread(fid, nt * nx, 'double=>double'), nt, nx);
  fclose(fid);
  if ~exist(output_dir, 'dir')
    mkdir(output_dir);
  end

  for k = 1:warmups
    run_psm(data, params);
  end
  samples = zeros(1, repeats);
  result = [];
  for k = 1:repeats
    started = tic;
    result = run_psm(data, params);
    samples(k) = toc(started);
  end

  fid = fopen(fullfile(output_dir, 'matlab_samples.csv'), 'wt');
  fprintf(fid, '%.9f', samples(1));
  fprintf(fid, ',%.9f', samples(2:end));
  fprintf(fid, '\n');
  fclose(fid);

  values = reshape(result.', [], 1);
  raw = zeros(2 * numel(values), 1);
  raw(1:2:end) = real(values);
  raw(2:2:end) = imag(values);
  fid = fopen(fullfile(output_dir, 'matlab_result.c64'), 'wb', 'ieee-le');
  fwrite(fid, raw, 'double');
  fclose(fid);
  fprintf('MATLAB PSM: %d samples\n', repeats);
end

function image = run_psm(data, p)
  nt = p.dims(1);
  nx = p.dims(2);
  fs = p.fs(1);
  delay = p.tDelay(1);
  velocity = p.cc(1);
  f_low = p.fLow(1);
  f_high = p.fHigh(1);
  x_step = p.xStep(1);
  thickness = p.thick(1);
  nfft_t = 2 ^ nextpow2(nt);
  nfft_x = 2 ^ nextpow2(nx);

  signed_t = 0:(nfft_t - 1);
  signed_t(signed_t >= ceil(nfft_t / 2)) = signed_t(signed_t >= ceil(nfft_t / 2)) - nfft_t;
  omega = signed_t * (2 * pi * fs / nfft_t);
  keep = omega >= 2 * pi * f_low & omega <= 2 * pi * f_high;
  omega = omega(keep);
  spectrum = fft(data, nfft_t, 1);
  spectrum = spectrum(keep, :) .* exp(-1i * omega * delay).';
  spectrum = fft(spectrum, nfft_x, 2).';

  signed_x = 0:(nfft_x - 1);
  signed_x(signed_x >= ceil(nfft_x / 2)) = signed_x(signed_x >= ceil(nfft_x / 2)) - nfft_x;
  kx = signed_x.' * ((2 * pi / x_step) / nfft_x);
  kz2 = (2 / velocity) ^ 2 * (ones(nfft_x, 1) * (omega .^ 2)) - (kx .^ 2) * ones(1, numel(omega));
  propagating = kz2 >= 0;
  kz = sqrt(max(kz2, 0));
  spectrum(~propagating) = 0;

  dz = (velocity / 2) / f_high;
  z_offset = delay * (velocity / 2);
  nz = ceil((thickness - z_offset) / dz) + 1;
  wave = spectrum .* exp(1i * kz * z_offset);
  phase = exp(1i * kz * dz);
  image = complex(zeros(nz, nx));
  for iz = 1:nz
    plane = ifft(sum(wave, 2));
    image(iz, :) = plane(1:nx).';
    wave = wave .* phase;
  end
end

function params = read_params(path)
  params = struct();
  lines = regexp(strtrim(fileread(path)), '\r?\n', 'split');
  for k = 1:numel(lines)
    parts = strsplit(strtrim(lines{k}));
    if numel(parts) > 1
      params.(parts{1}) = sscanf(strjoin(parts(2:end), ' '), '%f').';
    end
  end
end