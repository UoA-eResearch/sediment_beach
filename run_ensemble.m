function run_ensemble(params, S0, out_root)
% RUN_ENSEMBLE  Run ShorelineS for several parameter sets in parallel (parfor).
%
%   run_ensemble(params, S0, out_root)
%
%   params    struct array; each element holds the S fields to override, e.g.
%                 params = struct('qscal', {0.5, 1, 2});
%                 params = struct('qscal', {0.5, 1}, 'wvcfile', {'sinthe_waves.nc', 'input_25_contour.nc'});
%   S0        base settings (as in hindcast_run.m, with S0.xmc / S0.ymc set)
%   out_root  folder for the runs; run k writes to out_root/run_00k
%
% Each run writes output.mat and params.mat to its own folder. With the
% Parallel Computing Toolbox, parfor runs the parameter sets on the workers of
% the parallel pool (parpool(n) sets their number); without it, and in GNU
% Octave, the loop runs serially. Plots and videos are switched off.
% See python/run_ensemble.py for the Python version.

if nargin < 3 || isempty(out_root)
    out_root = 'ensemble';
end
if ~exist(out_root, 'dir')
    mkdir(out_root);
end
addpath(fullfile(fileparts(mfilename('fullpath')), 'ShorelineS_functions'));

% hidden figures (Octave's fltk toolkit cannot print a hidden figure, so it keeps them visible)
S0.plotvisible = double(exist('OCTAVE_VERSION','builtin') > 0);
S0.video = 0;
S0.plotinterval = 1e9;   % no plotting during the run
names = fieldnames(params);
n = numel(params);
status = cell(n,1);
wall = NaN(n,1);

parfor k = 1:n
    S = S0;
    for j = 1:numel(names)
        S.(names{j}) = params(k).(names{j});
    end
    S.outputdir = fullfile(out_root, sprintf('run_%03d', k));
    t0 = tic;
    try
        ShorelineS(S);
        status{k} = 'ok';
    catch err
        status{k} = ['failed: ' err.message];
    end
    wall(k) = toc(t0);
    save_params(fullfile(S.outputdir, 'params.mat'), params(k));
end

% summary
fid = fopen(fullfile(out_root, 'summary.csv'), 'w');
fprintf(fid, 'run,%s,wall_s,status\n', strjoin(names', ','));
for k = 1:n
    vals = cellfun(@(f) value_str(params(k).(f)), names', 'UniformOutput', false);
    fprintf(fid, '%d,%s,%.1f,%s\n', k, strjoin(vals, ','), wall(k), strrep(status{k}, ',', ';'));
end
fclose(fid);
fprintf('%d runs done (%d failed); summary: %s\n', n, sum(~strcmp(status,'ok')), ...
    fullfile(out_root, 'summary.csv'));
end


function save_params(fname, p)
% (a separate function: save cannot be called directly inside parfor)
if ~exist(fileparts(fname), 'dir')
    mkdir(fileparts(fname));
end
save(fname, '-struct', 'p');
end


function s = value_str(v)
if ischar(v)
    s = v;
elseif isnumeric(v) || islogical(v)
    s = strtrim(sprintf('%g ', v));
else
    s = class(v);
end
end
