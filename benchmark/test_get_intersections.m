function test_get_intersections(origdir)
% TEST_GET_INTERSECTIONS  Check that get_intersections gives bit-identical results to the original
%   origdir : folder with the original (unoptimised) ShorelineS functions, e.g. a
%             checkout of the commit before the optimisation:
%             git worktree add ../baseline aa45896  ->  origdir = '../baseline/ShorelineS_functions'
repo = fileparts(fileparts(mfilename('fullpath')));
O = origdir; N = fullfile(repo,'ShorelineS_functions');
try, rng(2); catch, rand('seed',2); randn('seed',2); end
nbad=0; tO=0; tN=0; ncr=0;
for trial=1:3000
  m=randi(60); n=randi(60);
  xi=cumsum(randn(1,m+1)); yi=cumsum(randn(1,m+1)); xj=cumsum(randn(1,n+1))+randn*3; yj=cumsum(randn(1,n+1));
  if rand<.3, xi(randi(m+1))=NaN; xj(randi(n+1))=NaN; end
  if rand<.3, k=randi(m); xi(k+1)=xi(k); end   % vertical segment
  if rand<.3, k=randi(n); xj(k+1)=xj(k); end
  if rand<.2, xi=round(xi); yi=round(yi); xj=round(xj); yj=round(yj); end % touching/colinear
  self = rand<.2;
  addpath(O); tic; if self, [a{1:8}]=get_intersections(xi,yi); else [a{1:8}]=get_intersections(xi,yi,xj,yj); end; tO=tO+toc; rmpath(O);
  addpath(N); tic; if self, [b{1:8}]=get_intersections(xi,yi); else [b{1:8}]=get_intersections(xi,yi,xj,yj); end; tN=tN+toc; rmpath(N);
  ncr=ncr+numel(a{1});
  for k=1:8
    if ~isequal(size(a{k}),size(b{k})) || ~isequaln(a{k},b{k})
      nbad=nbad+1; fprintf('MISMATCH trial %d out %d sizes %s %s\n',trial,k,mat2str(size(a{k})),mat2str(size(b{k}))); break;
    end
  end
end
% long coastline vs short ray (typical model use)
xc=cumsum(rand(1,700))*75; yc=200*sin(xc/3000);
tic; addpath(O); for r=1:300, [p{1:8}]=get_intersections(xc,yc,[5000+r 5000+r],[-500 500]); end; t1=toc; rmpath(O);
tic; addpath(N); for r=1:300, [q{1:8}]=get_intersections(xc,yc,[5000+r 5000+r],[-500 500]); end; t2=toc; rmpath(N);
fprintf('mismatches=%d crossings=%d random: old %.2fs new %.2fs | coast-vs-ray: old %.3fs new %.3fs (%.0fx) same=%d\n',nbad,ncr,tO,tN,t1,t2,t1/t2,isequal(p,q));
end
