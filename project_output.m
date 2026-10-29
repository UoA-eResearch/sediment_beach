%function project_output(x, y, file_out, dy, interv, x_initial, y_initial, n_d, wave_mat)
% PROJECT_OUTPUT  Compare the modelled and the observed (CoastSat/PCA) shorelines.
%
%   project_output(x, y, file_out, dy, interv, x_initial, y_initial, n_d)
%   project_output(..., wave_mat)
%
% Both coastlines are projected onto the initial coastline (normals of
% +-Lcrit) with get_polydistance and plotted as cross-shore change.
%
%   x, y        observed coastlines from initial_grid (time x transects,
%               daily, row 1 = 1 August 1999, the first CoastSat date)
%   file_out    folder with the model's output.mat (relative to pwd, or absolute)
%   dy          days between 1 August 1999 and the start of the simulation
%   interv      storage interval of the model output [days]
%   n_d         number of days simulated
%   wave_mat    wave data for the wave figures (default 'wave_data.mat', from
%               create_nc.m; use 'waves_sinth.mat' for the synthetic waves of
% %               sinth_waves.m). The wave figures are skipped if it cannot be read.
clear all
[x, y, x_initial, y_initial, g] = initial_grid;
dy = 0;
interv=30; %in days
wave_mat = 'wave_data_av.mat';
n_d = 9284;
% Dates are datenums so that this also runs in GNU Octave (no datetime)
date_start = datenum(1999,8,1);         % first CoastSat date (row 1 of x, y)

% Model output
output_file = fullfile('C:\Users\mpul348\Documents\sediment_beach\python\output\28_oct_per_per_tanbeta0033', 'output.mat');
if ~exist(output_file, 'file')
    output_file = fullfile(pwd, file_out, 'output.mat');
end
raw_out = load(output_file);
x_model = raw_out.O.x;   % transects x time
y_model = raw_out.O.y;

% Observed shoreline data over the simulated period
% (row dy+1 is the start date; dy = 0 when the run starts on 1 August 1999)
rows = (dy+1):min(dy+n_d+1, size(x,1));
x_obs = x(rows,:);   % time x transects
y_obs = y(rows,:);

% select observations every 'interv' days
if interv > 1
    idx_obs = 1:interv:size(x_obs,1);
    x_obst = x_obs(idx_obs,:)';
    y_obst = y_obs(idx_obs,:)';
else
    x_obst = x_obs';
    y_obst = y_obs';
end

%% Apply the function get_polydistance
Lcrit = 500;

% Number of points along the projection/grid
nGrid = size(x_obs,2);

% Number of time steps available in BOTH datasets
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

for it = 1:ntime

    % MODEL
    [dmin,xcr,ycr] = get_polydistance( ...
        x_initial,y_initial, ...
        x_model(:,it),y_model(:,it), ...
        Lcrit);

    zg_model(:,it) = dmin(:);
    xc_model(:,it) = xcr(:);
    yc_model(:,it) = ycr(:);

    % OBSERVED (interval)
    [dmin,xcr,ycr] = get_polydistance( ...
        x_initial,y_initial, ...
        x_obst(:,it),y_obst(:,it), ...
        Lcrit);

    zg_obs(:,it) = dmin(:);
    xc_obs(:,it) = xcr(:);
    yc_obs(:,it) = ycr(:);

end

%% Daily observed distance (independent of 'interv')

for it = 1:ntime_d

    [dmin,xcr,ycr] = get_polydistance( ...
        x_initial,y_initial, ...
        x_obst_d(:,it),y_obst_d(:,it), ...
        Lcrit);

    zg_obs_d(:,it) = dmin(:);
    xc_obs_d(:,it) = xcr(:);
    yc_obs_d(:,it) = ycr(:);

end

%% Plot


% Dates of the stored model output / interval observations
dates_out = date_start + interv*(0:ntime-1);
dates_d   = date_start + (0:ntime_d-1);
tr = 1:nGrid;   % transect index (alongshore position)

figure
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

nexttile
imagesc(dates_out, tr, zg_model)
set(gca,'YDir','normal');
datetick('x','mmm/yyyy','keeplimits')
title('Model')
xlabel('Time')
ylabel('Alongshore position (transect)')
caxis([-200 200])
colorbar

nexttile
imagesc(dates_out, tr, zg_obs)
set(gca,'YDir','normal');
datetick('x','mmm/yyyy','keeplimits')
title('CoastSat')
xlabel('Time')
caxis([-20 20])
colorbar

nexttile
sm=movmean(zg_obs_d,720,2,'omitnan');
imagesc(dates_out, tr, sm)
set(gca,'YDir','normal');
datetick('x','mmm/yyyy','keeplimits')
title('CoastSat (720-day moving mean)')
xlabel('Time')
caxis([-20 20])
colorbar
colormap(flipud(cmocean('balance')))

%% Observed and modelled shorelines, one per year (every 12th stored output)
explained_file = 'explained_rem_site2.mat';
variance_pct = NaN;

if exist(explained_file,'file')
    ev = load(explained_file, 'explained_rem');
    variance_pct = ev.explained_rem{3} * 100;
end

idx = 1:48:ntime;
dv = datevec(dates_out(idx));
years_plot = dv(:,1);
cmap = cmocean('haline', max(length(idx),2));
% 
% figure('Color','w','Position',[100 100 1000 900])
% tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
% 
% % ---- Observed shorelines ----
% ax1 = nexttile;
% hold(ax1,'on')
% for j = 1:length(idx)
%     plot(ax1, x_obst(:,idx(j)), y_obst(:,idx(j)), ...
%         'Color', cmap(j,:), 'LineWidth', 3);
% end
% title(ax1, sprintf('Observed shorelines (%.1f%% variance explained)', variance_pct))
% ylabel(ax1,'Y (m)')
% axis(ax1,'equal')
% box(ax1,'on')
% grid(ax1,'on')
% % Dummy scatter maps year values to the colormap for the colorbar
% scatter(ax1, nan, nan, 1, nan, 'filled');
% colormap(ax1, cmap);
% clim(ax1, [min(years_plot) max(years_plot)]);
% cb1 = colorbar(ax1);
% cb1.Label.String = 'Year';
% cb1.Ticks = years_plot(1: max(1, round(length(years_plot)/6)):end);
% cb1.TickLabels = string(cb1.Ticks);
% xlim(ax1,[1.894470e6 1.900640e6])
% 
% % ---- Modelled shorelines ----
% ax2 = nexttile;
% hold(ax2,'on')
% for j = 1:length(idx)
%     plot(ax2, x_model(:,idx(j)), y_model(:,idx(j)), ...
%         'Color', cmap(j,:), 'LineWidth', 3);
% end
% title(ax2,'Modelled shorelines')
% xlabel(ax2,'X (m)')
% ylabel(ax2,'Y (m)')
% axis(ax2,'equal')
% box(ax2,'on')
% grid(ax2,'on')
% scatter(ax2, nan, nan, 1, nan, 'filled');
% colormap(ax2, cmap);
% clim(ax2, [min(years_plot) max(years_plot)]);
% cb2 = colorbar(ax2);
% cb2.Label.String = 'Year';
% cb2.Ticks = years_plot(1: max(1, round(length(years_plot)/6)):end);
% cb2.TickLabels = string(cb2.Ticks);
% xlim(ax2,[1.894470e6 1.900640e6])
% 

%% Wave figures

a = load(wave_mat);
t  = a.data.dates;t = datenum(t);
ll = a.data.ysave;
h  = a.data.hs_int_av;     % time x nodes
w  = a.data.dpm_int_av;

% only the simulated period
date_end = dates_out(end);
in = t >= date_start & t <= date_end;
t = t(in); h = h(in,:); w = w(in,:);

% Alongshore component of the wave energy flux relative to the mean direction
d = mod(atan2d(mean(sind(w),1,'omitnan'), mean(cosd(w),1,'omitnan')), 360);  % 1 x nodes
sxy = h.^2 .* cosd(w-d) .* sind(w-d);

% figure('Color','w', 'Units','centimeters', 'Position',[2 2 12 24]);
% tiledlayout(4,1,'TileSpacing','compact','Padding','compact')
% 
% ax1 = nexttile;
% sxym=movmean(sxy,1080,1,'omitnan');
% pcolor(ax1,t,ll,sxym');shading flat;
% colormap(ax1,cmocean('curl'))
% desv_pad = std(sxym(~isnan(sxym)));
% caxis(ax1, [-2 2] * desv_pad);
% cb1 = colorbar(ax1);
% cb1.Ticks = [-1 0 1] * desv_pad;
% cb1.TickLabelInterpreter = 'latex';
% cb1.TickLabels = {'$-1\sigma$', '$0$', '$+1\sigma$'};
% ylabel(cb1,'S_{xy}')
% xlim(ax1,[date_start date_end])
% % set(ax1,'XTickLabel',[])
% ylabel(ax1,'Latitude (wave nodes)');
% xlim(ax1,[date_start date_end])
% datetick(ax1,'x','mmm/yyyy','keeplimits')
% 
% ax2 = nexttile;
% pcolor(ax2,dates_out,tr,zg_obs);
% shading flat;
% colormap(ax2,flipud(cmocean('balance')))
% caxis(ax2,[-20 20])
% cb2 = colorbar(ax2);
% ylabel(cb2,'CoastSat (m)')
% xlim(ax2,[date_start date_end])
% set(ax2,'XTickLabel',[])
% ylabel(ax2,'Transect');set(ax2, 'Color', 'k');  % NaN values appear black
% 
% ax4 = nexttile;
% sm=movmean(zg_obs_d,720,2,'omitnan');
% imagesc(ax4,dates_out,tr,sm);
% set(ax4, 'YDir', 'normal')
% colormap(ax4,flipud(cmocean('balance')))
% cb4 = colorbar(ax4);
% ylabel(cb4,'CoastSat (720-day moving mean)')
% xlim(ax4,[date_start date_end])
% datetick(ax4,'x','mmm/yyyy','keeplimits')
% ylabel(ax4,'Transect');
% set(ax4, 'Color', 'k');  % NaN values appear black
% caxis(ax4,[-10 10])
% 
% ax3 = nexttile;
% pcolor(ax3,dates_out,tr,zg_model);
% shading flat;
% colormap(ax3,flipud(cmocean('balance')))
% cb3 = colorbar(ax3);
% ylabel(cb3,'Model (m)')
% xlim(ax3,[date_start date_end])
% datetick(ax3,'x','mmm/yyyy','keeplimits')
% ylabel(ax3,'Transect');
% set(ax3, 'Color', 'k');  % NaN values appear black
% caxis(ax3,[-200 200])




%%%%%%%%%%
%%% RMSE 

% Define the transects you want to analyse
tr_select = 1:450;      

% Find their row indices in tr
idx_tr = ismember(tr, tr_select);

% Subset transects
tr_sub = tr(idx_tr);
zg_obs_sub = zg_obs(idx_tr,:);
zg_model_sub = zg_model(idx_tr,:);
sm_sub = sm(idx_tr,:);

%% Difference between observed and modelled shoreline change
err = zg_obs_sub - zg_model_sub;  % Selected transects x time
MAE = mean(abs(err(:,:)), 'all', 'omitnan');
rmse_transect = sqrt(mean(err.^2, 2, 'omitnan'));
rmse_time = sqrt(mean(err.^2, 1, 'omitnan'));
MAE = mean(abs(err), 'all', 'omitnan');
mae_transect = mean(abs(err), 2, 'omitnan');

figure('Color','w')
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

% 1. Shoreline error
ax1 = nexttile;
imagesc(ax1, dates_out, tr_sub, err);
xlim(ax1,[date_start date_end]);
datetick(ax1,'x','mmm/yyyy','keeplimits');
set(ax1,'YDir','normal');
colormap(ax1, flipud(cmocean('curl')));
cb1 = colorbar(ax1);
ylabel(cb1,'Observed - Model (m)');
ylabel(ax1,'Transect');
title(ax1, sprintf('Difference between observed and modelled'));
text(ax1, 0.02, 0.95, sprintf('MAE = %.2f m', MAE), ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'BackgroundColor','w', ...
    'EdgeColor','k', ...
    'Margin',4);
%caxis(ax1,[-80 80])

ax3 = nexttile;
imagesc(ax3, dates_out, tr_sub, sm_sub);
xlim(ax3,[date_start date_end]);
datetick(ax3,'x','mmm/yyyy','keeplimits');
set(ax3,'YDir','normal');
colormap(ax3, flipud(cmocean('balance')));
cb1 = colorbar(ax3);
caxis(ax3,[-10 10]);

ax4 = nexttile;
imagesc(ax4, dates_out, tr_sub, zg_model_sub);
xlim(ax4,[date_start date_end]);
datetick(ax4,'x','mmm/yyyy','keeplimits');
set(ax4,'YDir','normal');
colormap(ax4, flipud(cmocean('balance')));
cb1 = colorbar(ax4);
%caxis(ax4,[-40 40]);

% 2. RMSE by transect
ax2 = nexttile;
plot(ax2, rmse_transect, tr_sub, 'k-', 'LineWidth',1.2);
xlabel(ax2,'rmse (m)');
ylabel(ax2,'Transect');
title(ax2,'RMSE by transect');
set(ax2,'YDir','normal');
grid(ax2,'on');
