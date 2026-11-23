% Single-run wall-clock timings of the ORIGINAL Octave code on the tasks timed in
% examples/figures.rs (indicative only; one run each).
addpath('../synaptus/matlab/core');
ds = '../synaptus/datasets/'; fid = fopen('figures/data/timings_octave.tsv','w');
s = load([ds 'LineScan2D_WireTargets.mat']); [nT,~]=size(s.ptx); th=(s.tDelay+nT/s.fs)*s.cc/2;
tic; psm(s.ptx,s.fs,s.tDelay,s.cc,th,0.4e6,2.5e6,s.xStep); fprintf(fid,'wire PSM (psm.m)\t%.4f\n',toc);
s = load([ds 'PlaneScan3D_PlexiAluFBH.mat']);
tic; psm(s.ptxy,s.fs,s.tDelay,s.cc,s.thick,s.fLow,s.fHigh,[s.xStep s.yStep]); fprintf(fid,'plane 3D PSM (psm.m)\t%.4f\n',toc);
s = load([ds 'LineScan2D_PinsPlexiAluSDH.mat']);
tic; psm(s.ptx,s.fs,s.tDelay,s.cc,s.thick,s.fLow,s.fHigh,s.xStep); fprintf(fid,'pins PSM (psm.m)\t%.4f\n',toc);
tic; mulok(s.ptx,s.fs,s.tDelay,s.cc,s.thick,s.fLow,s.fHigh,s.xStep); fprintf(fid,'pins MULOK (mulok.m)\t%.4f\n',toc);
load([ds 'CylScan3D_ObjectsInPipe.mat'],'ptpz','fs','phiStep','zStep','tDelay','r0','fLow','fHigh','cc');
rRaw = r0 + (tDelay+(1:size(ptpz,1))/fs)*(cc/2); ptpz = ptpz(rRaw<=0.147,:,:);
rStep=(cc/2)/(fHigh-fLow); rStart=r0+tDelay*cc/2; rEnd=rStart+(size(ptpz,1)/fs)*(cc/2);
tic; cpsm(ptpz,fs,tDelay,cc,fLow,fHigh,[phiStep zStep],r0,rStep,rStart,rEnd,'transFunc','haun'); fprintf(fid,'pipe CPSM haun (cpsm.m)\t%.4f\n',toc);
fclose(fid);
