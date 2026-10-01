function t = bench_run(codedir, ndays, doplot, outdir)
% BENCH_RUN  Time a short ShorelineS hindcast with the settings of hindcast_run.m
%
%   t = bench_run(codedir, ndays, doplot, outdir)
%
%   codedir : folder with the ShorelineS functions to time
%             (e.g. '../ShorelineS_functions', or a copy of an older version)
%   ndays   : number of simulated days, starting 2000-01-01 (dt = 3 h)
%   doplot  : 1 = plot every time step (ShorelineS default plotinterval=1),
%             0 = plot only the first time step
%   outdir  : output folder for output.mat
%
% The initial coastline is taken from the stored run in 30_sept_int_w/output.mat,
% so initial_grid.m (which needs readtable and the Mapping Toolbox) is not needed.
% Works in MATLAB and in GNU Octave (with the netcdf package).

repo = fileparts(fileparts(mfilename('fullpath')));
addpath(codedir);
if exist('OCTAVE_VERSION','builtin')
    pkg load netcdf
end
r = load(fullfile(repo,'30_sept_int_w','output.mat'),'S');

S = struct;
S.reftime = '2000-01-01';
S.endofsimulation = datestr(datenum(2000,1,1)+ndays,'yyyy-mm-dd');
S.xmc = r.S.xmc;
S.ymc = r.S.ymc;
S.wvcfile = fullfile(repo,'input_25_contour.nc');
S.interpolationmethod = 'alongshoremapping';
S.ddeep = 25;
S.dnearshore = 25;
S.ds0 = 75;
S.trform = 'KAMP';
S.boundaryconditionstart = 'periodic';
S.boundaryconditionend = 'periodic';
S.d = 10;
S.dt = 1/(365*8);
S.tc = 0;
S.storageinterval = 30;
S.outputdir = outdir;
S.suppresshighangle = 0;
S.plotvisible = 1;
S.plotDIR = 30;
S.plotHS = 30;
S.video = 0;
S.fignryear = 1;
S.fastplot = 0;
if ~doplot
    S.plotinterval = 1e9;   % only the first time step is plotted
end
if ~exist(outdir,'dir')
    mkdir(outdir);
end

tic;
S = ShorelineS(S); %#ok<NASGU>
t = toc;
fprintf('BENCH total %.2f s\n', t);
rmpath(codedir);
end
