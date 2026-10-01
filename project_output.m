function project_output(x, y,file_out,dy,interv, x_initial, y_initial, n_d)

% Model output
output_folder = fullfile(pwd, file_out);
raw_out = load(fullfile(output_folder, 'output.mat'));

x_model = raw_out.O.x;   % alongshore × time
y_model = raw_out.O.y;

% Observed shoreline data
x_obs = x(dy:dy+n_d,:);   % time × transects
y_obs = y(dy:dy+n_d,:);

%select observations every 'interv' days

if interv > 1
    idx_obs = 1:interv:size(x_obs,1);
    x_obst = x_obs(idx_obs,:)';
    y_obst = y_obs(idx_obs,:)';
else
    x_obst = x_obs';
    y_obst = y_obs';
end

%%Apply the function get_polydistance
Lcrit = 500;

%Number of points along the projection/grid
nGrid = size(x_obs,2);

%Number of time steps available in BOTH datasets
ntime = min(size(x_model,2), size(x_obst,2));

% Preallocate
zg_model = NaN(nGrid,ntime);
xc_model = NaN(nGrid,ntime);
yc_model = NaN(nGrid,ntime);

zg_obs = NaN(nGrid,ntime);
xc_obs = NaN(nGrid,ntime);
yc_obs = NaN(nGrid,ntime);

%% Daily observed data
% Always calculate this using ALL observations, regardless of 'interv'

x_obst_d = x_obs';
y_obst_d = y_obs';

ntime_d = size(x_obst_d,2);

zg_obs_d = NaN(nGrid,ntime_d);
xc_obs_d = NaN(nGrid,ntime_d);
yc_obs_d = NaN(nGrid,ntime_d);

%% Calculate distances

for t = 1:ntime

    % ==========================
    % MODEL
    % ==========================
    [dmin,xcr,ycr] = get_polydistance( ...
        x_initial,y_initial, ...
        x_model(:,t),y_model(:,t), ...
        Lcrit);

    zg_model(:,t) = dmin(:);
    xc_model(:,t) = xcr(:);
    yc_model(:,t) = ycr(:);


    % ==========================
    % OBSERVED (interval)
    % ==========================
    [dmin,xcr,ycr] = get_polydistance( ...
        x_initial,y_initial, ...
        x_obst(:,t),y_obst(:,t), ...
        Lcrit);

    zg_obs(:,t) = dmin(:);
    xc_obs(:,t) = xcr(:);
    yc_obs(:,t) = ycr(:);

end

%% Daily observed distance
% This is independent of 'interv'

for t = 1:ntime_d

    [dmin,xcr,ycr] = get_polydistance( ...
        x_initial,y_initial, ...
        x_obst_d(:,t),y_obst_d(:,t), ...
        Lcrit);

    zg_obs_d(:,t) = dmin(:);
    xc_obs_d(:,t) = xcr(:);
    yc_obs_d(:,t) = ycr(:);

end

%% Plot
figure
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile
imagesc(zg_model)
set(gca,'YDir','reverse')
colorbar
title('Model')
xlabel('Time (yr)')
ylabel('Alongshore position')

nexttile
imagesc(1:size(zg_model,2), 1:size(zg_model,1), zg_obs_d)
set(gca,'YDir','reverse')
colorbar
title('CoastSat')
xlabel('Time (yr)')

nexttile
window = round(30); 
sm = NaN(size(zg_obs_d));
for i = 1:size(zg_obs_d,1)
    sm(i,:) = movmean(zg_obs_d(i,:), window, 'omitnan');
end
imagesc(1:size(zg_model,2), 1:size(zg_model,1), sm)
set(gca,'YDir','reverse')
colorbar
title('CoastSat (variability >2 yr )')
xlabel('Time (yr)')

colormap(cmocean('haline'))