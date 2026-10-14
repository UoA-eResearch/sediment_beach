clear all

% This function creates a NetCDF file with synthetic waves
% with synthetic wave directions

sinth = 'sinthe_waves.nc';

% Reading the NetCDF file from Oceanum database
ncfile = 'oceanum_paper2.nc';

lat      = ncread(ncfile,'latitude');  lat = lat(:);
lon      = ncread(ncfile,'longitude'); lon = lon(:);
time_raw = ncread(ncfile,'time');
hs       = ncread(ncfile,'hs');
dpm      = ncread(ncfile,'dpm');
depth    = ncread(ncfile,'botl');
tps      = ncread(ncfile,'tps');

% Visualize Oceanum nodes location
figure
scatter(lat,lon,30,mean(depth,2),'filled')
colorbar
xlabel('Latitude')
ylabel('Longitude')

% Select dates to match shoreline data
dy = 0;

time_origin = datetime(1979,2,1,0,0,0);
time_datetime = time_origin + hours(time_raw);

date_request = datetime(1999,8,1,0,0,0);
start_date_request = date_request + days(dy);

end_date_request = datetime(2024,12,30,0,0,0);

actual_start = min(time_datetime);
actual_end   = max(time_datetime);

start_date = max(start_date_request,actual_start);
end_date   = min(end_date_request,actual_end);

time_mask = time_datetime >= start_date & ...
            time_datetime <= end_date;

selected_idx = find(time_mask);

dates = time_datetime(selected_idx);
dates = dates(:);

hs_idx  = hs(:,selected_idx);
dpm_idx = dpm(:,selected_idx);
tps_idx = tps(:,selected_idx);

%% Find the 25 m contour

llat = min(lat):0.01:max(lat);
llon = min(lon):0.01:max(lon);

[LAT,LON] = meshgrid(llat,llon);

DEP = griddata(lat,lon,mean(depth,2),LAT,LON);

figure
pcolor(llon,llat,DEP')
colorbar
shading flat
xlabel('Longitude')
ylabel('Latitude')
colormap(flipud(cmocean('haline')))

M = contourc(llat,llon,DEP,[20 25 30]);

[j,k] = find(M(1,:) == 25);

xsave = NaN;
ysave = NaN;

for j = 1:length(k)

    hold on

    x = M(1,k(j)+1:k(j)+M(2,k(j)));
    y = M(2,k(j)+1:k(j)+M(2,k(j)));

    plot(y,x,'k')

    xsave = [xsave;x'];
    ysave = [ysave;y'];

end

xsave = xsave(2:end);
ysave = ysave(2:end);

%% Writing time vector in seconds

t_ref = datetime(1999,8,1,0,0,0);

time_sec = seconds(dates - t_ref);
Nt = length(time_sec);

nccreate(sinth,'time', ...
    'Dimensions',{'time',Nt}, ...
    'Datatype','double');

ncwrite(sinth,'time',time_sec);

%% Writing coordinates of wave nodes

Ns = length(ysave);

crs = projcrs(2193);

[station_x,station_y] = projfwd(crs,xsave,ysave);

nccreate(sinth,'station_x', ...
    'Dimensions',{'station',Ns}, ...
    'Datatype','double');

nccreate(sinth,'station_y', ...
    'Dimensions',{'station',Ns}, ...
    'Datatype','double');

ncwrite(sinth,'station_x',station_x);
ncwrite(sinth,'station_y',station_y);

%% Interpolation weights

Nn = numel(lon);

F_w = scatteredInterpolant( ...
    lon,lat,zeros(Nn,1), ...
    'linear','nearest');

W = zeros(Ns,Nn);

for k = 1:Nn

    e_k = zeros(Nn,1);
    e_k(k) = 1;

    F_w.Values = e_k;

    W(:,k) = F_w(ysave,xsave);

end

%% Interpolate Hs along 25 m contour

hs_int = W*hs_idx;

for it = find(any(isnan(hs_idx),1))

    F_w.Values = hs_idx(:,it);

    hs_int(:,it) = F_w(ysave,xsave);

end

%% Interpolate wave period along 25 m contour

tps_int = W*tps_idx;

for it = find(any(isnan(tps_idx),1))

    F_w.Values = tps_idx(:,it);

    tps_int(:,it) = F_w(ysave,xsave);

end

%% Create synthetic wave directions

rng(1)

n_years = floor(days(dates(end)-dates(1))/365.25) + 1;

dpm_int = NaN(Ns,Nt);

for iy = 1:n_years

    start_date = dates(1) + calyears(iy-1);

    if iy < n_years
        end_date = dates(1) + calyears(iy);
    else
        end_date = dates(end) + seconds(1);
    end

    idx_year = dates >= start_date & dates < end_date;

    n = sum(idx_year);

    if mod(iy,2) == 1

        % Odd years: directions between 0-45 or 225-360

        r = rand(Ns,n);

        idx1 = r < 0.5;

        tmp = NaN(Ns,n);

        tmp(idx1) = 45 * rand(sum(idx1(:)),1);

        tmp(~idx1) = 225 + 135 * rand(sum(~idx1(:)),1);

        dpm_int(:,idx_year) = tmp;

    else

        % Even years: directions between 45-225

        dpm_int(:,idx_year) = ...
            45 + 180 * rand(Ns,n);

    end

end

%% Convert to time x station

hs_int_m  = hs_int';
dpm_int_m = dpm_int';
tps_int_m = tps_int';

%% Write variables to NetCDF

nccreate(sinth,'point_hm0', ...
    'Dimensions',{'time',Nt,'station',Ns}, ...
    'Datatype','double');

ncwrite(sinth,'point_hm0',hs_int_m);

nccreate(sinth,'point_wavdir', ...
    'Dimensions',{'time',Nt,'station',Ns}, ...
    'Datatype','double');

ncwrite(sinth,'point_wavdir',dpm_int_m);

nccreate(sinth,'point_tp', ...
    'Dimensions',{'time',Nt,'station',Ns}, ...
    'Datatype','double');

ncwrite(sinth,'point_tp',tps_int_m);

%% NetCDF attributes

ncwriteatt(sinth,'time', ...
    'units','seconds since 1999-08-01 00:00:00');

ncwriteatt(sinth,'time', ...
    'long_name','time');

ncwriteatt(sinth,'point_hm0', ...
    'units','m');

ncwriteatt(sinth,'point_tp', ...
    'units','s');

ncwriteatt(sinth,'point_wavdir', ...
    'units','degrees');

%% Save variables

data.dates      = dates;
data.xsave      = xsave;
data.ysave      = ysave;
data.dpm_int_m  = dpm_int_m;
data.hs_int_m   = hs_int_m;
data.tps_int_m  = tps_int_m;

save('waves_sinth.mat','data','-v7.3')

%% Plot results

a = load('waves_sinth.mat');

t  = a.data.dates;
ll = a.data.ysave;
h  = a.data.hs_int_m;
w  = a.data.dpm_int_m;
p  = a.data.tps_int_m;

figure

subplot(3,1,1)

pcolor(t,ll,w')
shading flat
set(gca,'YDir','reverse')

ylabel('Longitude')

ax = gca;
ax.XTickLabel = [];

colormap(flipud(cmocean('phase',256)))
clim([0 360])

cb2 = colorbar;
cb2.Ticks = [0 90 180 270 360];
cb2.TickLabels = {'0','90','180','270','360'};
ylabel(cb2,'Wave direction (°)')


subplot(3,1,2)

pcolor(t,ll,h')
shading flat
set(gca,'YDir','reverse')

ylabel('Longitude')

ax = gca;
ax.XTickLabel = [];

cmocean('haline')

cb1 = colorbar;
ylabel(cb1,'H_s (m)')


subplot(3,1,3)

pcolor(t,ll,p')
shading flat
set(gca,'YDir','reverse')

ylabel('Longitude')

ax = gca;
ax.XAxis.TickLabelFormat = 'yyyy';
ax.XTickLabel = string(year(ax.XTick));

cmocean('thermal')

cb3 = colorbar;
ylabel(cb3,'T_p (s)')

xlabel('Timestep = 3 hours')