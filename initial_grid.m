

function [x, y, x_initial, y_initial, g] = initial_grid(beach, dy, n_d)
% ==============================================================
% INITIAL GRID
%
% Reconstructs CoastSat shoreline positions using the PCA model
% and restores transects that were removed from the PCA.
%
% Outputs:
%   x, y          = time × transects, with NaN separators
%   x_initial     = initial shoreline x coordinates
%   y_initial     = initial shoreline y coordinates
%   g             = [x_initial y_initial]
%
% Important:
%   x and y contain the FULL time series.
%   dy and n_d are not used to trim the output here.
%   Trimming is done later in project_output.
% ==============================================================


%% --------------------------------------------------------------
% File names
% ---------------------------------------------------------------

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


%% --------------------------------------------------------------
% Transects removed from the PCA reconstruction
% ---------------------------------------------------------------

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


%% --------------------------------------------------------------
% Load PCA shoreline data
% ---------------------------------------------------------------

A_file = 'A_site2.mat';
mint_file = 'mint_site2.mat';

load(A_file,'A');
load(mint_file,'mint');


%% --------------------------------------------------------------
% Coordinate reference system
% ---------------------------------------------------------------

crs = projcrs(2193);       % NZTM2000


%% --------------------------------------------------------------
% Initialise combined arrays
% ---------------------------------------------------------------

x = [];
y = [];

x_initial = [];
y_initial = [];


%% ==============================================================
% Loop through selected beaches
% ==============================================================

for ib = 1:length(beach)

    b = beach(ib);

    fprintf('\nProcessing beach %d...\n',b);


    %% ----------------------------------------------------------
    % Read transect information
    % -----------------------------------------------------------

    T = readtable(fileNames{b},'Sheet','Transects');

    % Longitude
    land_x = T{:,18};

    % Latitude
    land_y = T{:,19};

    % Transect orientation
    phi = T{:,3};


    %% ----------------------------------------------------------
    % Number of transects
    % -----------------------------------------------------------

    nTransects = length(land_x);


    %% ----------------------------------------------------------
    % Transects retained for PCA reconstruction
    % -----------------------------------------------------------

    idx_remove = remove{b};

    idx_valid = true(nTransects,1);

    idx_valid(idx_remove) = false;


    %% ----------------------------------------------------------
    % Convert transect origins to NZTM2000
    % -----------------------------------------------------------

    [land_x_utm,land_y_utm] = projfwd( ...
        crs,land_y,land_x);


    %% ----------------------------------------------------------
    % Transect direction vectors
    % -----------------------------------------------------------

    theta = deg2rad(phi);

    ux = sin(theta);
    uy = cos(theta);


    %% ----------------------------------------------------------
    % Coordinates and directions for retained transects
    % -----------------------------------------------------------

    P0 = [land_x_utm(idx_valid), ...
        land_y_utm(idx_valid)];


    ux_valid = ux(idx_valid);
    uy_valid = uy(idx_valid);


    %% ----------------------------------------------------------
    % PCA shoreline reconstruction
    % -----------------------------------------------------------

    A_s = A{b};

    mint_s = mint{b};

    % Convert PCA anomaly to absolute shoreline position
    A_abs = A_s + mint_s;


    %% ----------------------------------------------------------
    % Check dimensions
    % -----------------------------------------------------------

    if size(A_abs,2) ~= sum(idx_valid)

        error(['Beach %d: number of PCA transects (%d) ' ...
               'does not match number of valid transects (%d).'], ...
               b,size(A_abs,2),sum(idx_valid));

    end


    %% ----------------------------------------------------------
    % Reconstruct shoreline coordinates
    % -----------------------------------------------------------

    x_valid = ...
        P0(:,1)' + A_abs .* ux_valid';

    y_valid = ...
        P0(:,2)' + A_abs .* uy_valid';


    % Number of time steps
    nTime = size(x_valid,1);


    %% ==========================================================
    % Read satellite observations for removed transects
    % ===========================================================

    time_obs = readtable( ...
        fileNames{b}, ...
        'Sheet','Intersect points', ...
        'Range','A2:A517');

    s_obs = readtable( ...
        fileNames{b}, ...
        'Sheet','Intersect points', ...
        'Range','B2:KH517');


    %% ----------------------------------------------------------
    % Convert table to cell array
    % -----------------------------------------------------------

    data = table2cell(s_obs);

    [nRows,nCols] = size(data);


    %% ----------------------------------------------------------
    % Prepare latitude and longitude arrays
    % -----------------------------------------------------------

    lat = NaN(nRows,nCols);
    lon = NaN(nRows,nCols);


    %% ----------------------------------------------------------
    % Parse "lat lon" strings
    % -----------------------------------------------------------

    for i = 1:nRows

        for j = 1:nCols

            val = data{i,j};

            % Empty cell
            if isempty(val)
                continue
            end

            % Numeric NaN
            if isnumeric(val) && isnan(val)
                continue
            end

            % Convert to string
            val = string(val);

            % Replace spaces with commas
            val = strrep(val," ",",");

            % Split
            parts = split(val,",");

            % Need at least latitude and longitude
            if numel(parts) >= 2

                lat(i,j) = str2double(parts(1));

                lon(i,j) = str2double(parts(2));

            end

        end

    end


    %% ----------------------------------------------------------
    % Convert satellite coordinates to NZTM2000
    % -----------------------------------------------------------

    [x_obs,y_obs] = projfwd(crs,lat,lon);


    %% ==========================================================
    % Create complete shoreline matrix for this beach
    % ===========================================================

    x_b = NaN(nTime,nTransects);

    y_b = NaN(nTime,nTransects);


    %% ----------------------------------------------------------
    % Insert PCA-reconstructed transects
    % -----------------------------------------------------------

    x_b(:,idx_valid) = x_valid;

    y_b(:,idx_valid) = y_valid;


    %% ==========================================================
    % Restore removed transects
    % ==========================================================
    %
    % For each removed transect:
    %
    %   1. Find the first valid satellite observation
    %   2. Use that coordinate
    %   3. Keep that coordinate constant through time
    %
    % This means the removed transects are NOT reconstructed
    % using PCA.
    % ===========================================================

    for tr = idx_remove(:)'

        % Column in Intersect points corresponding to this
        % original transect
        obs_col = tr;


        % Check that the column exists
        if obs_col > size(x_obs,2)

            warning(['Beach %d, transect %d: no corresponding ' ...
                     'satellite column found.'],b,tr);

            continue

        end


        % Find first valid satellite position
        valid_obs = ...
            ~isnan(x_obs(:,obs_col)) & ...
            ~isnan(y_obs(:,obs_col));


        if any(valid_obs)

            first_idx = find(valid_obs,1,'first');

            x_removed = x_obs(first_idx,obs_col);
            y_removed = y_obs(first_idx,obs_col);


            % Keep the coordinate constant through all model time
            x_b(:,tr) = x_removed;
            y_b(:,tr) = y_removed;

        else

            warning(['Beach %d, transect %d: no valid satellite ' ...
                     'observation found.'],b,tr);

        end

    end


    %% ----------------------------------------------------------
    % Initial shoreline for this beach
    % -----------------------------------------------------------

    x_initial_b = x_b(1,:);

    y_initial_b = y_b(1,:);


    %% ==========================================================
    % Combine beaches
    % ===========================================================

    if isempty(x)

        % -------------------------------------------------------
        % First beach
        % -------------------------------------------------------

        x = x_b;

        y = y_b;

        x_initial = x_initial_b;

        y_initial = y_initial_b;


    else

        % -------------------------------------------------------
        % Subsequent beaches
        %
        % Add ONE NaN column between beaches
        % -------------------------------------------------------

        x = [ ...
            x, ...
            NaN(size(x,1),1), ...
            x_b];


        y = [ ...
            y, ...
            NaN(size(y,1),1), ...
            y_b];


        x_initial = [ ...
            x_initial, ...
            NaN, ...
            x_initial_b];


        y_initial = [ ...
            y_initial, ...
            NaN, ...
            y_initial_b];

    end


    fprintf('  Total transects: %d\n',nTransects);

    fprintf('  PCA transects:   %d\n',sum(idx_valid));

    fprintf('  Removed:         %d\n',length(idx_remove));

end


%% ==============================================================
% Final initial coastline
% ==============================================================

g = [ ...
    x_initial(:), ...
    y_initial(:)];


%% --------------------------------------------------------------
% Diagnostics
% ---------------------------------------------------------------

fprintf('\n========================================\n');
fprintf('INITIAL GRID COMPLETE\n');
fprintf('========================================\n');

fprintf('Selected beaches: ');

fprintf('%d ',beach);

fprintf('\n');

fprintf('x dimensions:          %d × %d\n',size(x,1),size(x,2));

fprintf('y dimensions:          %d × %d\n',size(y,1),size(y,2));

fprintf('x_initial dimensions:  %d × %d\n',size(x_initial,1),size(x_initial,2));

fprintf('Number of separators:  %d\n',sum(isnan(x_initial)));

fprintf('Separator locations:   ');

disp(find(isnan(x_initial)));

fprintf('========================================\n');

