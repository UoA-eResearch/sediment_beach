function test_get_polydistance(origdir)
% TEST_GET_POLYDISTANCE  Compare vectorised get_polydistance with the original loop (real + random data)
%   origdir : folder with the original (unoptimised) ShorelineS functions, e.g. a
%             checkout of the commit before the optimisation:
%             git worktree add ../baseline aa45896  ->  origdir = '../baseline/ShorelineS_functions'
repo = fileparts(fileparts(mfilename('fullpath')));
O = origdir; N = fullfile(repo,'ShorelineS_functions');
try, rng(3); catch, rand('seed',3); randn('seed',3); end
r=load(fullfile(repo,'30_sept_int_w','output.mat'));
Xr=r.S.xmc; Yr=r.S.ymc; worst=0; nm=0; tO=0; tN=0;
cols=[1 10 25 52];
for c=cols
  Xc=r.O.x(:,c); Yc=r.O.y(:,c);
  addpath(O); tic; [a1,a2,a3,a4,a5,a6]=get_polydistance(Xr,Yr,Xc,Yc,500); tO=tO+toc; rmpath(O);
  addpath(N); tic; [b1,b2,b3,b4,b5,b6]=get_polydistance(Xr,Yr,Xc,Yc,500); tN=tN+toc; rmpath(N);
  A={a1,a2,a3,a4,a5,a6}; B={b1,b2,b3,b4,b5,b6};
  for k=1:6
    if ~isequal(isnan(A{k}),isnan(B{k})), nm=nm+1; fprintf('NaN pattern differs col %d out %d (%d vs %d)\n',c,k,sum(isnan(A{k})),sum(isnan(B{k}))); end
    d=max(abs(A{k}(:)-B{k}(:))); if ~isempty(d), worst=max(worst,d); end
  end
  fprintf('col %d: valid=%d identical=%d\n',c,sum(~isnan(a1)),isequaln(A,B));
end
% random polylines incl. NaN separators
for t=1:300
  n=randi(80)+2; Xr=cumsum(rand(1,n))*50; Yr=cumsum(randn(1,n))*10;
  m=randi(120)+2; Xc=cumsum(rand(1,m))*40-100; Yc=cumsum(randn(1,m))*10+randn*100;
  if rand<.3, Xc(randi(m))=NaN; Yc(Xc~=Xc)=NaN; end
  if rand<.3, Xr(randi(n))=NaN; Yr(Xr~=Xr)=NaN; end
  addpath(O); [a1,a2,a3,a4,a5,a6]=get_polydistance(Xr,Yr,Xc,Yc,200); rmpath(O);
  addpath(N); [b1,b2,b3,b4,b5,b6]=get_polydistance(Xr,Yr,Xc,Yc,200); rmpath(N);
  A={a1,a2,a3,a4,a5,a6}; B={b1,b2,b3,b4,b5,b6};
  for k=1:6
    if ~isequal(isnan(A{k}),isnan(B{k})), nm=nm+1; fprintf('rand %d out %d NaN pattern differs\n',t,k); end
    d=max(abs(A{k}(:)-B{k}(:))); if ~isempty(d), worst=max(worst,d); end
  end
end
fprintf('NaN mismatches=%d worst absdiff=%g | real data: old %.2fs new %.3fs (%.0fx)\n',nm,worst,tO,tN,tO/tN);
end
