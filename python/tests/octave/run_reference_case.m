function run_reference_case(infile, outfile, fundir)
% RUN_REFERENCE_CASE  Run ShorelineS (functions in fundir) with the settings S in
% infile and save its output structure O to outfile (-v7), for the Python tests.
% Called by make_references.py.
addpath(fundir);
pkg load netcdf
s = load(infile);
S = s.S;
if isfield(S, 'boundaryconditionstart') && iscell(S.boundaryconditionstart)
    S.boundaryconditionstart = reshape(S.boundaryconditionstart, 1, []);
end
if isfield(S, 'boundaryconditionend') && iscell(S.boundaryconditionend)
    S.boundaryconditionend = reshape(S.boundaryconditionend, 1, []);
end
S.xmc = S.xmc(:)';
S.ymc = S.ymc(:)';
% no plots or videos
S.plotvisible = 1; S.video = 0; S.fignryear = 1; S.fastplot = 0; S.plotinterval = 1e9;
if ~exist(S.outputdir, 'dir'), mkdir(S.outputdir); end
ShorelineS(S);
r = load(fullfile(S.outputdir, 'shorelines_output.mat'));
O = r.O;
save('-v7', outfile, 'O');
end
