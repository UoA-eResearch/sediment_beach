function make_function_references(datadir, fundir)
% MAKE_FUNCTION_REFERENCES  Results of single ShorelineS functions (in fundir) for
% the cases in datadir, for the Python tests. Called by make_references.py.
addpath(fundir);
warning('off', 'all');

% get_intersections: 400 random cases
c = load(fullfile(datadir, 'get_intersections_cases.mat')); C = c.C;
R = cell(size(C,1), 8);
for k = 1:size(C,1)
    if C{k,5}
        [R{k,1:8}] = get_intersections(C{k,1}, C{k,2});
    else
        [R{k,1:8}] = get_intersections(C{k,1}, C{k,2}, C{k,3}, C{k,4});
    end
end
save('-v7', fullfile(datadir, 'get_intersections_octave.mat'), 'R');

% get_clockpoly: 300 random polygons
c = load(fullfile(datadir, 'clockpoly_cases.mat')); C = c.C;
R = zeros(size(C,1), 1);
for k = 1:size(C,1)
    R(k) = get_clockpoly(C{k,1}, C{k,2});
end
save('-v7', fullfile(datadir, 'clockpoly_octave.mat'), 'R');

% merge_coastlines_mc: 400 random pairs of crossing sections
c = load(fullfile(datadir, 'merge_crossing_input.mat')); cases = c.cases;
out = cell(size(cases,1), 2);
for k = 1:size(cases,1)
    COAST = struct('x_mc', cases{k,1}(:)', 'y_mc', cases{k,2}(:)', 'ds0', 75);
    evalc('COAST = merge_coastlines_mc(COAST);');
    out(k,:) = {COAST.x_mc, COAST.y_mc};
end
save('-v7', fullfile(datadir, 'merge_crossing_octave.mat'), 'out');

% merge_coastlines_mc: the merge of two beaches in the hindcast with synthetic waves
c = load(fullfile(datadir, 'merge_event_synthetic_waves.mat')); ev = c.ev;
out = cell(size(ev,1), 2);
for k = 1:size(ev,1)
    COAST = struct('x_mc', ev{k,1}(:)', 'y_mc', ev{k,2}(:)', 'ds0', 75);
    evalc('COAST = merge_coastlines_mc(COAST);');
    out(k,:) = {COAST.x_mc, COAST.y_mc};
end
save('-v7', fullfile(datadir, 'merge_event_synthetic_waves_octave.mat'), 'out');

% find_overwash_mc: a thin closed island whose first point is shadowed
d = load(fullfile(datadir, 'overwash_cyclic_case.mat'));
COAST = struct('x_mc', d.x_mc(:)', 'y_mc', d.y_mc(:)', 'n_mc', 1, 'ds0', 75, ...
               'PHIcxy_mc', d.PHIcxy_mc(:)');
nx = numel(COAST.x_mc); nq = numel(d.QS_mc);
COAST.s_mc = zeros(1, nx); COAST.ds_mc = zeros(1, nx); COAST.h0_mc = 10*ones(1, nx);
COAST.PHIcxy0_mc = COAST.PHIcxy_mc;
WAVE = struct('HSo_mc', ones(1, nq), 'TP_mc', ones(1, nq), 'PHItdp_mc', d.PHIbr_mc(:)', ...
              'PHIbr_mc', d.PHIbr_mc(:)', 'diff_mc', false(1, nq));
TRANSP = struct('trform', 'KAMP', 'QS_mc', d.QS_mc(:)', 'shadow_mc', logical(d.shadow_mc(:)'), ...
                'shadowS_hD_mc', false(1, nq));
SPIT = struct('method', 'default', 'spitwidth', 50, 'owscale', 0.1, 'owtimescale', 0, 'Dsf', 8, ...
              'Dbb', 4, 'bheight', 2);
STRUC = struct('xhard', [], 'yhard', []);
DUNE = struct('used', 0); MUD = struct('used', 0); TIME = struct('dt', 1/(365*8));
COAST = find_overwash_mc(COAST, WAVE, SPIT, STRUC, TRANSP, DUNE, MUD, TIME);
x_mc = COAST.x_mc; y_mc = COAST.y_mc;
save('-v7', fullfile(datadir, 'overwash_cyclic_octave.mat'), 'x_mc', 'y_mc');
end
