% Export Synaptus datasets to raw little-endian binary and produce reference
% outputs from the ORIGINAL psm.m / mulok.m (legacy behaviour) for the Rust
% port's consistency check. Run from synaptus-rs/:  octave --no-gui octave/export_and_reference.m
% Output format: <name>.f64 = raw float64 LE, MATLAB column-major order;
%                <name>.txt = "key value..." lines (dims, params).
% Complex arrays are written interleaved (re, im) per element.
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
  for i = 1:numel(f)
    v = s.(f{i}); fprintf(fid, '%s', f{i}); fprintf(fid, ' %.17g', double(v(:))); fprintf(fid, '\n');
  end
  fclose(fid);
end
function export_result(prefix, im, zIm)
  for l = 1:numel(im)
    write_cplx(sprintf('%s_layer%d.c64', prefix, l), im{l});
    s = struct('dims', size(im{l}), 'z', zIm{l}); write_params(sprintf('%s_layer%d.txt', prefix, l), s);
  end
end

addpath('../synaptus/matlab/core');
ds = '../synaptus/datasets/'; out = 'data/';

% ---- 2D single layer: WireTargets (band/thickness as in learn_psm.m / Python test)
s = load([ds 'LineScan2D_WireTargets.mat']);
[nT, nX] = size(s.ptx); fLow = 0.4e6; fHigh = 2.5e6;
thick = (s.tDelay + nT/s.fs) * s.cc/2;
write_real([out 'wire.f64'], s.ptx);
write_params([out 'wire.txt'], struct('dims',[nT nX 1],'fs',s.fs,'tDelay',s.tDelay,'cc',s.cc, ...
  'thick',thick,'fLow',fLow,'fHigh',fHigh,'xStep',s.xStep,'yStep',s.xStep));
[im,~,zIm] = psm(s.ptx, s.fs, s.tDelay, s.cc, thick, fLow, fHigh, s.xStep); export_result([out 'wire_psm'], im, zIm);
[im,~,zIm] = mulok(s.ptx, s.fs, s.tDelay, s.cc, thick, fLow, fHigh, s.xStep); export_result([out 'wire_mulok'], im, zIm);

% ---- 2D three layers: PinsPlexiAluSDH
s = load([ds 'LineScan2D_PinsPlexiAluSDH.mat']); [nT, nX] = size(s.ptx);
write_real([out 'pins.f64'], s.ptx);
write_params([out 'pins.txt'], struct('dims',[nT nX 1],'fs',s.fs,'tDelay',s.tDelay,'cc',s.cc, ...
  'thick',s.thick,'fLow',s.fLow,'fHigh',s.fHigh,'xStep',s.xStep,'yStep',s.xStep));
[im,~,zIm] = psm(s.ptx, s.fs, s.tDelay, s.cc, s.thick, s.fLow, s.fHigh, s.xStep); export_result([out 'pins_psm'], im, zIm);
[im,~,zIm] = mulok(s.ptx, s.fs, s.tDelay, s.cc, s.thick, s.fLow, s.fHigh, s.xStep); export_result([out 'pins_mulok'], im, zIm);

% ---- 3D three layers: PlaneScan3D_PlexiAluFBH
s = load([ds 'PlaneScan3D_PlexiAluFBH.mat']); [nT, nX, nY] = size(s.ptxy);
write_real([out 'plane.f64'], s.ptxy);
write_params([out 'plane.txt'], struct('dims',[nT nX nY],'fs',s.fs,'tDelay',s.tDelay,'cc',s.cc, ...
  'thick',s.thick,'fLow',s.fLow,'fHigh',s.fHigh,'xStep',s.xStep,'yStep',s.yStep));
[im,~,~,zIm] = psm(s.ptxy, s.fs, s.tDelay, s.cc, s.thick, s.fLow, s.fHigh, [s.xStep s.yStep]); export_result([out 'plane_psm'], im, zIm);
[im,~,~,zIm] = mulok(s.ptxy, s.fs, s.tDelay, s.cc, s.thick, s.fLow, s.fHigh, [s.xStep s.yStep]); export_result([out 'plane_mulok'], im, zIm);
printf('\nOK\n');
