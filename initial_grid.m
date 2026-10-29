


function [x, y, x_initial, y_initial, g] = initial_grid
% ==============================================================
% INITIAL GRID
%
% Reconstruct shoreline positions using the PCA model.
% Transects excluded from the PCA are completely removed from the
% returned x, y, x_initial, y_initial and g arrays.
% A single NaN separator is retained between different beaches.
%
% Outputs:
%   x, y       = time x retained transects, with NaN separators
%                between beaches
%   x_initial  = initial shoreline x coordinates, with separators
%   y_initial  = initial shoreline y coordinates, with separators
%   g          = [x_initial(:), y_initial(:)]
%
% dy and n_d are retained for compatibility with the original
% function signature. Trimming is handled later in project_output.
% ==============================================================
beach = [2 3];

fileNames = { ...
    'nzd0204.xlsx'
    'nzd0207.xlsx'
    'nzd0217.xlsx'
    'nzd0222.xlsx'
    'nzd0220.xlsx'
    'nzd0226.xlsx'
    'nzd0229.xlsx'
    'nzd0231.xlsx'
    'nzd0234.xlsx'
    'nzd0239.xlsx'
    'nzd0240.xlsx'
    'nzd0238.xlsx'
    'nzd0236.xlsx'
    'nzd0233.xlsx'
    'nzd0230.xlsx'
    'nzd0227.xlsx'};

%% Transects removed from PCA reconstruction
remove = cell(16,1);
remove{1}  = 92:112;
remove{2}  = [1:60, 233:293];
remove{3}  = 1:38;
remove{4}  = 36:43;
remove{5}  = 14:16;
remove{6}  = 165:166;
remove{7}  = 190:196;
remove{8}  = [1:5, 98:106];
remove{9}  = 145:150;
remove{10} = 110:141;
remove{11} = [1:136, 207:213, 254:264];
remove{12} = [];
remove{13} = 46:49;
remove{14} = [];
remove{15} = 26:60;
remove{16} = 27:43;

%% Load PCA shoreline data
A_file = 'PC_lst_site2.mat';
mint_file = 'mean_PC_lst_site2.mat';

load(A_file, 'PC_lst');
load(mint_file, 'mean_PC_lst');

%% Coordinate reference system
crs = projcrs(2193);  % NZTM2000

%% Initialise combined arrays
x = [];
y = [];
x_initial = [];
y_initial = [];

%% Process selected beaches
for ib = 1:numel(beach)

    b = beach(ib);
    fprintf('\nProcessing beach %d...\n', b);

    if b < 1 || b > numel(fileNames) || b ~= floor(b)
        error('Invalid beach index: %d.', b);
    end

    %% Read transect information
    T = readtable(fileNames{b}, 'Sheet', 'Transects');

    % Longitude, latitude and transect orientation
    land_x = T{:,18};
    land_y = T{:,19};
    phi = T{:,3};

    nTransects = numel(land_x);

    %% Identify retained transects
    idx_remove = remove{b};
    idx_remove = idx_remove(idx_remove >= 1 & idx_remove <= nTransects);

    idx_valid = true(nTransects,1);
    idx_valid(idx_remove) = false;

    %% Convert transect origins to NZTM2000
    [land_x_utm, land_y_utm] = projfwd(crs, land_y, land_x);

    %% Transect direction vectors
    theta = deg2rad(phi);
    ux = sin(theta);
    uy = cos(theta);

    %% Retained transect origins and direction vectors
    P0 = [land_x_utm(idx_valid), land_y_utm(idx_valid)];
    ux_valid = ux(idx_valid);
    uy_valid = uy(idx_valid);

    %% PCA shoreline reconstruction
    A_s = PC_lst{b};
    mint_s = mean_PC_lst{b};
    A_abs = A_s + mint_s;

    %% Check PCA dimensions
    if size(A_abs,2) ~= sum(idx_valid)
        error(['Beach %d: PCA has %d transect columns, but %d ' ...
            'transects are retained after exclusions.'], ...
            b, size(A_abs,2), sum(idx_valid));
    end

    %% Reconstruct retained shoreline coordinates
    x_b = P0(:,1)' + A_abs .* ux_valid';
    y_b = P0(:,2)' + A_abs .* uy_valid';


    %% Combine beaches, retaining one NaN separator between beaches
    if isempty(x)
        x = x_b;
        y = y_b;
        x_initial = x_b(1,:);
        y_initial = y_b(1,:);
    else
        x = [x, NaN(size(x,1),1), x_b];
        y = [y, NaN(size(y,1),1), y_b];

        x_initial = [x_initial, NaN, x_b(1,:)];
        y_initial = [y_initial, NaN, y_b(1,:)];
    end

    fprintf('  Original transects: %d\n', nTransects);
    fprintf('  PCA-retained:       %d\n', sum(idx_valid));
    fprintf('  Excluded/removed:   %d\n', numel(idx_remove));
end

%% Final initial coastline
% Retain separator rows (NaN, NaN) between beaches.
g = [x_initial(:), y_initial(:)];

%% Diagnostics
fprintf('\n========================================\n');
fprintf('INITIAL GRID COMPLETE\n');
fprintf('========================================\n');
fprintf('Selected beaches: ');
fprintf('%d ', beach);
fprintf('\n');
fprintf('x dimensions:          %d x %d\n', size(x,1), size(x,2));
fprintf('y dimensions:          %d x %d\n', size(y,1), size(y,2));
fprintf('x_initial dimensions:  %d x %d\n', size(x_initial,1), size(x_initial,2));
fprintf('Number of beach separators: %d\n', sum(isnan(x_initial)));
fprintf('========================================\n');

%% Plot initial coastline
figure('Color','w')
plot(g(:,1), g(:,2), 'k.-', ...
    'LineWidth', 0.8, ...
    'MarkerSize', 3);
axis equal;
box on;
xlabel('Easting (m)');
ylabel('Northing (m)');
title('Initial coastline — selected transects');
set(gca, 'FontSize', 12, 'TickDir', 'out');

