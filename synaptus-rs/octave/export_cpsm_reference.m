% Export ObjectsInPipe (cropped as in cpsm_test3D_ObjectsInPipe.m) and the
% ORIGINAL cpsm.m output ('haun' and 'gardner') for the T7 consistency check.
1;
function write_real(path, a)
  fid = fopen(path, 'w', 'ieee-le'); fwrite(fid, double(a(:)), 'float64'); fclose(fid);
end
function write_cplx(path, a)
  a = double(a(:)); v = zeros(2*numel(a),1); v(1:2:end) = real(a); v(2:2:end) = imag(a);
  fid = fopen(path, 'w', 'ieee-le'); fwrite(fid, v, 'float64'); fclose(fid);
end
function write_params(path, s)
  fid = fopen(path, 'w'); f = fieldnames(s);
  for i = 1:numel(f), v = s.(f{i}); fprintf(fid, '%s', f{i}); fprintf(fid, ' %.17g', double(v(:))); fprintf(fid, '\n'); end
  fclose(fid);
end
addpath('../synaptus/matlab/core');
load('../synaptus/datasets/CylScan3D_ObjectsInPipe.mat','ptpz','fs','phiStep','zStep','tDelay','r0','fLow','fHigh','cc');
rRaw = r0 + (tDelay+(1:size(ptpz,1))/fs) * (cc/2);
ptpz = ptpz(rRaw <= 0.147,:,:);
rStep = (cc/2)/(fHigh-fLow); rStart = r0 + (tDelay*(cc/2)); rEnd = rStart + (size(ptpz,1)/fs)*(cc/2);
rr = rStart:rStep:rEnd;
write_real('data/pipe.f64', ptpz);
write_params('data/pipe.txt', struct('dims', size(ptpz), 'fs', fs, 'tDelay', tDelay, 'fLow', fLow, 'fHigh', fHigh, ...
  'phiStep', phiStep, 'zStep', zStep, 'r0', r0, 'cc', cc, 'rr', rr));
for tf = {'haun','gardner'}
  [im,phiIm,rIm,zIm] = cpsm(ptpz,fs,tDelay,cc,fLow,fHigh,[phiStep zStep],r0,rStep,rStart,rEnd,'transFunc',tf{1});
  write_cplx(['data/pipe_cpsm_' tf{1} '.c64'], im);
  write_params(['data/pipe_cpsm_' tf{1} '.txt'], struct('dims', size(im)));
end
printf('OK\n');
