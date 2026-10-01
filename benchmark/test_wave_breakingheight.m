function test_wave_breakingheight(origdir)
% TEST_WAVE_BREAKINGHEIGHT  Compare vectorised wave_breakingheight with the original per-point loop
%   origdir : folder with the original (unoptimised) ShorelineS functions, e.g. a
%             checkout of the commit before the optimisation:
%             git worktree add ../baseline aa45896  ->  origdir = '../baseline/ShorelineS_functions'
repo = fileparts(fileparts(mfilename('fullpath')));
O = origdir; N = fullfile(repo,'ShorelineS_functions');
try, rng(1); catch, rand('seed',1); randn('seed',1); end
n=5000;
W.dPHItdp=[ (rand(1,n)-.5)*200, 89.9, -89.9, 0, 45, 90 ];
m=numel(W.dPHItdp);
W.HStdp=[rand(1,n)*6, 0.0001, 3, 0, 10, 1]; W.TP=[1+rand(1,n)*18, 2, 30, 8, 0.5, 12];
W.gamma=0.72; W.dnearshore=25; T.suppresshighangle=0; T.trform='KAMP';
addpath(O); tic; A=wave_breakingheight(W,T); t1=toc; rmpath(O);
addpath(N); tic; B=wave_breakingheight(W,T); t2=toc;
f={'HSbr','dPHIbr','hbr','cbr','nbr'};
for k=1:5, d=max(abs(A.(f{k})(:)-B.(f{k})(:))); fprintf('%s maxdiff=%g identical=%d size=%d\n',f{k},d,isequaln(A.(f{k}),B.(f{k})),isequal(size(A.(f{k})),size(B.(f{k})))); end
fprintf('old %.3fs new %.3fs speedup %.0fx\n',t1,t2,t1/t2);
rmpath(N);
end
