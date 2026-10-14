%function project_output(x, y,file_out,dy,interv, x_initial, y_initial, n_d)

% Model output
raw_out = load('C:\Users\mpul348\Documents\code_paper2_optimized\python\output\9_oct_closed\output.mat');
x_model = raw_out.O.x;   % transects × time
y_model = raw_out.O.y;   % transects × tim

% Observed shoreline data
x_obs = x; 
y_obs = y;

% %n_d is number of days in the simulation
% select observations every 'interv' days 

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
% 
% %% Plot

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
idx3=293+remove{3};

zg_model_plot=zg_model;
zg_obs_plot=zg_obs;

zg_model_plot([idx2 idx3],:)=NaN;
zg_obs_plot([idx2 idx3],:)=NaN;

sm=movmean(zg_obs_plot,30,2,'omitnan');
sm([idx2 idx3],:)=NaN;

figure
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')
n_t = (datetime(1999,8,1):calmonths(1):datetime(2004,6,1))';
n_t.Format = 'MMM/yyyy';
nexttile
imagesc(n_t, ll, zg_model_plot)
set(gca,'YDir','reverse','Color','w')
title('Model')
xlabel('Time')
ylabel('Alongshore position')
clim([-50 50])
colorbar

nexttile
imagesc(n_t, ll, zg_obs_plot)
set(gca,'YDir','reverse','Color','w')
title('CoastSat')
xlabel('Time')
colorbar
colormap(cmocean('haline'))


%%%
%plotting with waves

% a = load('wave_data.mat');
% t = a.data.dates;
% ll = a.data.ysave;
% la = a.data.xsave;
% h = a.data.hs_int_m;
% w = a.data.dpm_int_m;
% p = a.data.tps_int_m;
% 
% 
figure('Color','w')
tiledlayout(3,1,'TileSpacing','compact','Padding','compact')

ax1 = nexttile;
d = nan(size(w,2),1);
for si = 1:size(w,2)
    d(si) = atan2(mean(sin(w(:,si)/180*pi),'omitnan'), ...
                  mean(cos(w(:,si)/180*pi),'omitnan')) / pi * 180;
end
d = mod(d,360);
d = ones(size(w,1),1) * d';   % time × nodes
sxy = h.^2 .* cos((w-d)/180*pi) .* sin((w-d)/180*pi);%calculating the alongshore component relative to the mean direction

imagesc(ax1,t,ll,sxy')
set(ax1,'YDir','reverse','Color','w')
cmap = cmocean('balance',50);
colormap(ax1,cmap)
desv_pad = std(sxy(:),'omitnan');
lim = 2*desv_pad;
clim(ax1,[-lim lim])
cb1 = colorbar(ax1);
cb1.Ticks = [-1 0 1]*desv_pad;
cb1.TickLabels = {'-1\sigma','0','+1\sigma'};
cb1.Label.String = 'S_{xy}';
xticklabels(ax1,[])

ax1 = nexttile;
imagesc(t,ll,zg_obs_plot);
set(ax1,'YDir','reverse','Color','w')
cb2 = colorbar;
colormap(ax1,flipud(cmocean('curl')))
clim([-50 50])
xticklabels([])
ylabel('Longitude');
cb2.Label.String = 'CoastSat (m)';

ax2 = nexttile;
imagesc(t, ll,zg_model_plot);
set(ax1,'YDir','reverse','Color','w');
cb3 = colorbar;
cb3.Label.String = 'Model (m)';
colormap(ax2,flipud(cmocean('curl')))
%clim([-200 50])
xlabel('Time');

% 
% 

% 
% %%%%%
% %calculating annual cumulative sxy
% years = year(t);
% uyears = unique(years);
% sxy_year = NaN(length(uyears),size(sxy,2));
% for i = 1:length(uyears)
%     idx = years == uyears(i);
%     sxy_year(i,:) = sum(sxy(idx,:),1,'omitnan');
% end
% 
% cmap = cmocean('balance',56);
% lim = max(abs(sxy_year(:)),[],'omitnan');
% figure('Color','w')
% ax = axes;
% hold(ax,'on')
% dy_year = 0.003;
% for i = 1:length(uyears)
%     y_plot = y_wave + (i-1)*dy_year;
%     for j = 1:length(x_wave)-1
%         c = mean(sxy_year(i,j:j+1),'omitnan');
%         idx = round((c + lim)/(2*lim)*(size(cmap,1)-1)) + 1;
%         idx = max(1,min(size(cmap,1),idx));
%         plot(ax,x_wave(j:j+1),y_plot(j:j+1), ...
%             'Color',cmap(idx,:), ...
%            'LineWidth',3.2)
%     end
% end
% colormap(ax,cmap)
% xlabel(ax,'Longitude')
% ylabel(ax,'Latitude + year offset')
% cb = colorbar(ax);
% clim(ax,[-lim lim])
% cb.Label.String = 'Cumulative annual S_{xy}';
% outdir = 'C:\Users\mpul348\Documents\code_paper2_optimized\python\output\5_oct_opt';
% exportgraphics(gcf, ...
%     fullfile(outdir,'annual_cumulative_Sxy.png'), ...
%     'Resolution',1000);
% 

%ploting sinthetic waves along with modes

% a = load('waves_sinth.mat');
% t  = a.data.dates;
% ll = a.data.ysave;
% w  = a.data.dpm_int_m;
% start_date = datetime(1999,8,1);
% end_date   = datetime(2004,6,1);
% 
% idx = t >= start_date & t <= end_date;
% 
% t = t(idx);
% w = w(idx,:);
% 
% 
% figure('Color','w')
% tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
% 
% 
% %% 1. Wave direction
% 
% ax1 = nexttile;
% 
% imagesc(ax1,t,ll,w')
% 
% set(ax1,'YDir','reverse','Color','w')
% 
% colormap(ax1,cmocean('phase'))
% 
% clim(ax1,[0 360])
% 
% cb1 = colorbar(ax1);
% cb1.Label.String = 'Wave direction (°)';
% cb1.Ticks = [0 90 180 270 360];
% cb1.TickLabels = {'0','90','180','270','360'};
% 
% ylabel(ax1,'Longitude')
% xticklabels(ax1,[])
% 
% 
% %% 2. Model shoreline
% 
% ax2 = nexttile;
% 
% imagesc(ax2,t,ll,zg_model_plot)
% 
% set(ax2,'YDir','reverse','Color','w')
% 
% colormap(ax2,cmocean('haline'))
% 
% clim(ax2,[-50 50])
% 
% cb2 = colorbar(ax2);
% cb2.Label.String = 'Model (m)';
% 
% xlabel(ax2,'Time')
% 
% 
% xtickformat(ax2,'MMM/yyyy')




%% Plot observed and modelled shorelines
figure('Color','w','Position',[100 100 1000 900])
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
n_d = 64;
years_plot = 1999:2004;
idx = 1:12:n_d;  % August of each year
cmap = cmocean('haline',length(idx));
load('explained_rem_site2.mat', 'explained_rem');
variance_pct = explained_rem{3} * 100;

%% Observed shorelines
ax1 = nexttile;
hold(ax1,'on')

for j = 1:length(idx)
    i = idx(j);
    plot(ax1,x_obst(:,i),y_obst(:,i), ...
        'Color',cmap(j,:),'LineWidth',2);
end

title(ax1, sprintf('Observed shorelines (%.1f%% variance explained)', ...
    variance_pct))
ylabel(ax1,'Y (m)')
axis(ax1,'equal')
box(ax1,'on')
grid(ax1,'on')
colormap(ax1,cmap)
clim(ax1,[1 length(idx)])

cb1 = colorbar(ax1);
cb1.Label.String = 'Year';
cb1.Ticks = 1:length(idx);
cb1.TickLabels = string(years_plot);
xlim(ax1, [1.894e6 1.8951e6])

%% Modelled shorelines
ax1 = nexttile;
hold(ax1,'on')

for j = 1:length(idx)
    i = idx(j);
    plot(ax1,x_model(:,i),y_model(:,i), ...
        'Color',cmap(j,:),'LineWidth',2);
end
title(ax1,'Modelled shorelines')
xlabel(ax1,'X (m)')
ylabel(ax1,'Y (m)')
axis(ax1,'equal')
box(ax1,'on')
grid(ax1,'on')
colormap(ax1,cmap)
clim(ax1,[1 length(idx)])
cb2 = colorbar(ax1);
cb2.Label.String = 'Year';
cb2.Ticks = 1:length(idx);
cb2.TickLabels = string(years_plot);
xlim(ax1, [1.894e6 1.8951e6])
%%%%
%%%


%% Plot the PC WITHOUT the cross-displacement

a = load('wave_data.mat');
t = a.data.dates;
ll = a.data.ysave;
la = a.data.xsave;
h = a.data.hs_int_m;
w = a.data.dpm_int_m;
p = a.data.tps_int_m;

% figure('Color','w', 'Units','centimeters', ...
%        'Position',[2 2 12 20])
% tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
% load('explained_rem_site2.mat', 'explained_rem');
% variance_pct = explained_rem{3} * 100;
% 
% ax1 = nexttile;
% mean_dir = atan2d(mean(sind(w), 1, 'omitnan'), mean(cosd(w), 1, 'omitnan'));
% w_anom = mod(w - mean_dir + 180, 360) - 180;
% idx = t >= dates_plot(1) & t <= dates_plot(end);
% ll_we = flipud(ll);% Flip spatial order: west to east
% w_anom_we = fliplr(w_anom);
% pcolor(ax1, t(idx), ll_we, w_anom_we(idx,:)');
% shading(ax1, 'flat')
% colormap(ax1, cmocean('curl'))
% clim(ax1, [-50 50])
% xlim(ax1, [dates_plot(1) dates_plot(end)])
% set(ax1, 'YDir', 'normal')
% xticklabels([])
% cb1 = colorbar(ax1);
% cb1.Label.String = 'Wave direction anomaly (°)';
% cb1.Ticks = [-45 0 45];
% cb1.TickLabels = {'-45', '0', '45'};
% ylabel(ax1,'25m depth contour longitude (West to East)');
% 
% ax2 = nexttile;
% n_d = 132;
% dates_plot = datetime(1999,8,1) + calmonths(0:n_d-1);
% pcolor(ax2, 1:n_d, 1:size(zg_obs_plot,1), ...
%        zg_obs_plot(:,1:n_d))
% % Ticks every January, April, July, and October
% tick_months = [1 4 7 10];
% idx_ticks = ismember(month(dates_plot), tick_months);
% xticks(ax2, find(idx_ticks))
% xticklabels(ax2, cellstr(datestr(dates_plot(idx_ticks), 'mmm/yyyy')))
% xtickangle(ax2, 45)
% shading(ax2,'interp')
% colormap(ax2,flipud(cmocean('balance')))
% clim(ax2,[-20 20])
% cb2 = colorbar(ax2);
% ylabel(ax2, 'Transect index (West to East)')
% xlabel(ax2,'Date');
% xlim(ax2,[1 n_d])
% xtickangle(ax2,45)
% cb2.Label.String = sprintf('CoastSat (%.1f%% variance explained)', variance_pct);


% 
% Model
% ax2 = nexttile;
% pcolor(ax2, 1:n_d, 1:size(zg_model_plot,1), zg_model_plot(:,1:n_d))
% shading(ax2,'flat')
% colormap(ax2,flipud(cmocean('balance')))
% colorbar(ax2)
% ylabel(ax2,'Longitude');
% xlabel(ax2,'Date');
% title(ax2,'Model (m)')
% xlim(ax2,[1 n_d])
% xticks(ax2,1:12:n_d)
% xticklabels(ax2,datestr(dates_plot(1:12:n_d),'mmm-yyyy'))
% xtickangle(ax2,45)


% Choose plotting period
date_start = datetime(1999,8,1);
date_end   = datetime(2024,8,1);

figure('Color','w', 'Units','centimeters', ...
       'Position',[2 2 12 20])

tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

%% Panel 1: Wave direction anomalies
ax1 = nexttile;

mean_dir = atan2d(mean(sind(w), 1, 'omitnan'), ...
                  mean(cosd(w), 1, 'omitnan'));

w_anom = mod(w - mean_dir + 180, 360) - 180;

% Select dates
idx = t >= date_start & t <= date_end;

ll_we = flipud(ll);
w_anom_we = fliplr(w_anom);

pcolor(ax1, t(idx), ll_we, w_anom_we(idx,:)');
shading(ax1, 'flat')

colormap(ax1, cmocean('curl'))
clim(ax1, [-50 50])
xlim(ax1, [date_start date_end])

set(ax1, 'YDir', 'normal')
ylabel(ax1, '25 m depth contour longitude (West to East)')

% Quarterly ticks: Jan, Apr, Jul, Oct
tick_dates = (dateshift(date_start,'start','month'):calmonths(1): ...
              dateshift(date_end,'start','month'))';
tick_dates = tick_dates(ismember(month(tick_dates), [1 4 7 10]));

xticks(ax1, tick_dates)
xticklabels(ax1, cellstr(datestr(tick_dates,'mmm/yyyy')))
xtickangle(ax1, 45)

cb1 = colorbar(ax1);
cb1.Label.String = 'Wave direction anomaly (°)';
cb1.Ticks = [-45 0 45];

%% Panel 2: CoastSat observations
ax2 = nexttile;

n_d = size(zg_obs_plot, 2);
dates_coast = datetime(1999,8,1) + calmonths(0:n_d-1);

% Select dates
idx_coast = dates_coast >= date_start & dates_coast <= date_end;
dates_selected = dates_coast(idx_coast);

pcolor(ax2, dates_selected, 1:size(zg_obs_plot,1), ...
       zg_obs_plot(:,idx_coast));

shading(ax2, 'interp')
colormap(ax2, flipud(cmocean('balance')))
clim(ax2, [-20 20])

xlim(ax2, [date_start date_end])

% Use the same quarterly ticks
xticks(ax2, tick_dates)
xticklabels(ax2, cellstr(datestr(tick_dates,'mmm/yyyy')))
xtickangle(ax2, 45)

ylabel(ax2, 'Transect index (West to East)')
xlabel(ax2, 'Date')

cb2 = colorbar(ax2);
cb2.Label.String = sprintf( ...
    'CoastSat (%.1f%% variance explained)', variance_pct);


%%%%



% Choose plotting period
date_start = datetime(1999,8,1);
date_end   = datetime(2024,12,30);

% CoastSat explained variance
load('explained_rem_site2.mat', 'explained_rem');
variance_pct = explained_rem{3} * 100;

figure('Color','w', 'Units','centimeters', ...
       'Position',[2 2 12 20])

tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

%% Panel 1: Wave direction anomalies
ax1 = nexttile;

mean_dir = atan2d(mean(sind(w), 1, 'omitnan'), ...
                  mean(cosd(w), 1, 'omitnan'));

w_anom = mod(w - mean_dir + 180, 360) - 180;

% Select dates
idx = t >= date_start & t <= date_end;

% Spatial order: west to east
ll_we = flipud(ll);
w_anom_we = fliplr(w_anom);

pcolor(ax1, t(idx), ll_we, w_anom_we(idx,:)');
shading(ax1, 'flat')

colormap(ax1, cmocean('curl'))
clim(ax1, [-50 50])
xlim(ax1, [date_start date_end])

set(ax1, 'YDir', 'normal')
ylabel(ax1, '25 m depth contour longitude (West to East)')

cb1 = colorbar(ax1);
cb1.Label.String = 'Wave direction anomaly (°)';
cb1.Ticks = [-45 0 45];

%% Panel 2: CoastSat observations
ax2 = nexttile;

n_d = size(zg_obs_plot, 2);
dates_coast = datetime(1999,8,1) + calmonths(0:n_d-1);

% Select dates
idx_coast = dates_coast >= date_start & dates_coast <= date_end;
dates_selected = dates_coast(idx_coast);

pcolor(ax2, dates_selected, 1:size(zg_obs_plot,1), ...
       zg_obs_plot(:,idx_coast));

shading(ax2, 'interp')
colormap(ax2, flipud(cmocean('balance')))
clim(ax2, [-20 20])
xlim(ax2, [date_start date_end])

ylabel(ax2, 'Transect index (West to East)')
xlabel(ax2, 'Date')

cb2 = colorbar(ax2);
cb2.Label.String = sprintf( ...
    'CoastSat (%.1f%% variance explained)', variance_pct);

%% Yearly ticks: January, shared by both panels
tick_dates = (dateshift(date_start,'start','year'):calyears(1): ...
              dateshift(date_end,'start','year'))';

year_labels = cellstr(datestr(tick_dates,'yyyy'));

xticks(ax1, tick_dates)
xticklabels(ax1, year_labels)
xtickangle(ax1, 45)

xticks(ax2, tick_dates)
xticklabels(ax2, year_labels)
xtickangle(ax2, 45)


