function project_output(x, y, file_out, dy, interv, x_initial, y_initial, n_d, wave_mat)
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
%               sinth_waves.m). The wave figures are skipped if it cannot be read.

if nargin < 9 || isempty(wave_mat)
    wave_mat = 'wave_data.mat';
end
% Dates are datenums so that this also runs in GNU Octave (no datetime)
date0 = datenum(1999,8,1);         % first CoastSat date (row 1 of x, y)
date_start = date0 + dy;           % start of the simulation

% Model output
output_file = fullfile(file_out, 'output.mat');
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

% Transects removed from the PCA (held fixed in initial_grid) are not plotted.
% Beach 2 has 293 transects; beach 3 follows after one NaN separator.
remove = cell(16,1);
remove{1}=92:112;
remove{2}=[1:60,233:293];
remove{3}=1:38;
remove{4}=36:43;
remove{5}=14:16;
remove{6}=165:166;
remove{7}=190:196;
remove{8}=[1:5,98:106];
remove{9}=145:150;
remove{10}=110:141;
remove{11}=[1:136,207:213,254:264];
remove{12}=[];
remove{13}=46:49;
remove{14}=[];
remove{15}=26:60;
remove{16}=27:43;

idx2=remove{2};
idx3=293+1+remove{3};
idx_hide = [idx2 idx3];
idx_hide = idx_hide(idx_hide <= nGrid);

zg_model_plot=zg_model;
zg_obs_plot=zg_obs;

zg_model_plot(idx_hide,:)=NaN;
zg_obs_plot(idx_hide,:)=NaN;

sm=movmean(zg_obs_d,30,2,'omitnan');   % 30-day moving mean of the daily data
sm(idx_hide,:)=NaN;

% Dates of the stored model output / interval observations
dates_out = date_start + interv*(0:ntime-1);
dates_d   = date_start + (0:ntime_d-1);
tr = 1:nGrid;   % transect index (alongshore position)

figure
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

nexttile
imagesc(dates_out, tr, zg_model_plot)
set(gca,'YDir','reverse','Color','w')
datetick('x','mmm/yyyy','keeplimits')
title('Model')
xlabel('Time')
ylabel('Alongshore position (transect)')
caxis([-50 50])
colorbar

nexttile
imagesc(dates_out, tr, zg_obs_plot)
set(gca,'YDir','reverse','Color','w')
datetick('x','mmm/yyyy','keeplimits')
title('CoastSat')
xlabel('Time')
caxis([-50 50])
colorbar

nexttile
imagesc(dates_d, tr, sm)
set(gca,'YDir','reverse','Color','w')
datetick('x','mmm/yyyy','keeplimits')
title('CoastSat (30-day moving mean)')
xlabel('Time')
caxis([-50 50])
colorbar
colormap(cmocean('haline'))

%% Observed and modelled shorelines, one per year (every 12th stored output)
explained_file = 'explained_rem_site2.mat';
variance_pct = NaN;
if exist(explained_file,'file')
    ev = load(explained_file, 'explained_rem');
    variance_pct = ev.explained_rem{3} * 100;
end
idx = 1:12:ntime;
cmap = cmocean('haline',max(length(idx),2));
dv = datevec(dates_out(idx));
years_plot = dv(:,1);

figure('Color','w','Position',[100 100 1000 900])
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
hold(ax1,'on')
for j = 1:length(idx)
    plot(ax1,x_obst(:,idx(j)),y_obst(:,idx(j)),'Color',cmap(j,:),'LineWidth',2);
end
title(ax1, sprintf('Observed shorelines (%.1f%% variance explained)', variance_pct))
ylabel(ax1,'Y (m)')
axis(ax1,'equal'); box(ax1,'on'); grid(ax1,'on')
year_colorbar(ax1, cmap, years_plot)
xlim(ax1, [1.894e6 1.8951e6])

ax2 = nexttile;
hold(ax2,'on')
for j = 1:length(idx)
    plot(ax2,x_model(:,idx(j)),y_model(:,idx(j)),'Color',cmap(j,:),'LineWidth',2);
end
title(ax2,'Modelled shorelines')
xlabel(ax2,'X (m)')
ylabel(ax2,'Y (m)')
axis(ax2,'equal'); box(ax2,'on'); grid(ax2,'on')
year_colorbar(ax2, cmap, years_plot)
xlim(ax2, [1.894e6 1.8951e6])

%% Wave figures
try
    a = load(wave_mat);
    t  = a.data.dates;
    ll = a.data.ysave;
    h  = a.data.hs_int_m;     % time x nodes
    w  = a.data.dpm_int_m;
    if isa(t,'datetime')
        t = datenum(t);
    end
    t = t(:);
catch err
    warning('project_output:waves', ...
        'Wave figures skipped: cannot read %s (%s).', wave_mat, err.message);
    return
end

% only the simulated period
date_end = dates_out(end);
in = t >= date_start & t <= date_end;
t = t(in); h = h(in,:); w = w(in,:);

% Alongshore component of the wave energy flux relative to the mean direction
d = mod(atan2d(mean(sind(w),1,'omitnan'), mean(cosd(w),1,'omitnan')), 360);  % 1 x nodes
sxy = h.^2 .* cosd(w-d) .* sind(w-d);

figure('Color','w')
tiledlayout(3,1,'TileSpacing','compact','Padding','compact')

ax1 = nexttile;
imagesc(ax1,t,ll,sxy')
set(ax1,'YDir','reverse','Color','w')
colormap(ax1,cmocean('balance',50))
desv_pad = std(sxy(~isnan(sxy)));
caxis(ax1,[-2 2]*desv_pad)
cb1 = colorbar(ax1);
set_ticks(cb1,[-1 0 1]*desv_pad,{'-1\sigma','0','+1\sigma'});
ylabel(cb1,'S_{xy}')
xlim(ax1,[date_start date_end])
set(ax1,'XTickLabel',[])
ylabel(ax1,'Latitude (wave nodes)')

ax2 = nexttile;
imagesc(ax2,dates_out,tr,zg_obs_plot);
set(ax2,'YDir','reverse','Color','w')
colormap(ax2,flipud(cmocean('curl')))
caxis(ax2,[-50 50])
cb2 = colorbar(ax2);
ylabel(cb2,'CoastSat (m)')
xlim(ax2,[date_start date_end])
set(ax2,'XTickLabel',[])
ylabel(ax2,'Transect');

ax3 = nexttile;
imagesc(ax3,dates_out,tr,zg_model_plot);
set(ax3,'YDir','reverse','Color','w');
colormap(ax3,flipud(cmocean('curl')))
cb3 = colorbar(ax3);
ylabel(cb3,'Model (m)')
xlim(ax3,[date_start date_end])
datetick(ax3,'x','mmm/yyyy','keeplimits')
ylabel(ax3,'Transect');
xlabel(ax3,'Time');

% Wave direction anomalies and CoastSat change
mean_dir = atan2d(mean(sind(w),1,'omitnan'), mean(cosd(w),1,'omitnan'));
w_anom = mod(w - mean_dir + 180, 360) - 180;

figure('Color','w','Units','centimeters','Position',[2 2 12 20])
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

ax1 = nexttile;
pcolor(ax1, t, flipud(ll), fliplr(w_anom)');   % west to east
shading(ax1,'flat')
colormap(ax1, cmocean('curl'))
caxis(ax1,[-50 50])
xlim(ax1,[date_start date_end])
datetick(ax1,'x','yyyy','keeplimits')
set(ax1,'YDir','normal')
ylabel(ax1,'25 m depth contour (West to East)')
cb1 = colorbar(ax1);
ylabel(cb1,'Wave direction anomaly (°)')

ax2 = nexttile;
pcolor(ax2, dates_out, tr, zg_obs_plot);
shading(ax2,'interp')
colormap(ax2, flipud(cmocean('balance')))
caxis(ax2,[-20 20])
xlim(ax2,[date_start date_end])
datetick(ax2,'x','yyyy','keeplimits')
ylabel(ax2,'Transect index (West to East)')
xlabel(ax2,'Date')
cb2 = colorbar(ax2);
ylabel(cb2,sprintf('CoastSat (%.1f%% variance explained)', variance_pct))

end


function year_colorbar(ax, cmap, years_plot)
% Discrete colour bar labelled with the year of each plotted shoreline
n = size(cmap,1);
colormap(ax,cmap)
caxis(ax,[0.5 n+0.5])
cb = colorbar(ax);
ylabel(cb,'Year')
set_ticks(cb,1:numel(years_plot),cellstr(num2str(years_plot(:))));
end


function set_ticks(cb, ticks, labels)
% Colour bar ticks (MATLAB ColorBar properties; GNU Octave's colour bar is an axes)
if isprop(cb,'Ticks')
    set(cb,'Ticks',ticks,'TickLabels',labels);
else
    set(cb,'YTick',ticks,'YTickLabel',labels);
end
end
