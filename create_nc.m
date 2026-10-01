

clear all

%This function creates a netcdf file,input_nc, to be used in ShorelineS
%model by interpolating wave parametrs along a contour

input_nc = 'input_25_contour.nc'; %giving the file name

%reading the netcdf file from Oceanum database
ncfile = ('oceanum_paper2.nc'); 
lat  = ncread(ncfile, 'latitude'); lat = lat(:);
lon = ncread(ncfile, 'longitude');  lon = lon(:);
time_raw = ncread(ncfile, 'time');
hs  = ncread(ncfile, 'hs'); %height
dpm = ncread(ncfile, 'dpm'); %direction
depth = ncread(ncfile, 'botl');
tps = ncread(ncfile, 'tps'); %period

%visualize Oceanum nodes location
figure
scatter(lat,lon,30,mean(depth,2),'filled')
colorbar

%selecting the days on wave data to match the shoreline data 
dy=0; %days since the coastsat starts
time_origin = datetime(1979,2,1,0,0,0); %this is when the oceanum hindcast starts
time_datetime = time_origin + hours(time_raw); %time vector in hours
date_request = datetime(1999,08,01,0,0,0); 
start_date_request = date_request + days(dy);% shoreline data starts
end_date_request   = datetime(2024,12,30,0,0,0); %shoreline data ends
actual_start = min(time_datetime); actual_end = max(time_datetime); start_date = max(start_date_request, actual_start); end_date   = min(end_date_request, actual_end);
time_mask = (time_datetime >= start_date) & (time_datetime <= end_date);
selected_idx = find(time_mask); %index that matchs the shoreline data
dates = time_datetime(selected_idx); dates = dates(:);
hs_idx = hs(:, selected_idx);
dpm_idx = dpm(:, selected_idx);
tps_idx = tps(:, selected_idx);

% find the 20m contour
llat=min(lat):0.01:max(lat);
llon=min(lon):0.01:max(lon);
[LAT,LON]=meshgrid(llat,llon);
[m,n]=size(LAT); %getting grid dimensions
DEP=griddata(lat,lon,mean(depth,2),LAT,LON);
figure
pcolor(llon, llat, DEP');colorbar;shading flat;xlabel('Longitude');ylabel('Latitude');colormap(flipud(cmocean('haline')))
M=contourc(llat,llon,DEP,[20 25 30]);

[j,k]=find(M(1,:)==25)
xsave=NaN;
ysave=NaN;
for j=1:length(k) %"k" contains the positions where each 25 m contour starts. So we go through each 25 m contour one at a time.
    hold on
    x=M(1,k(j)+1:k(j)+M(2,k(j)));
    y=M(2,k(j)+1:k(j)+M(2,k(j)));
   plot(y,x,'k')

   xsave=[xsave;x'];
   ysave=[ysave;y'];
end
xsave=xsave(2:end);ysave=ysave(2:end) % those are the contour coordinates

%writing time vector in seconds

t_ref = datetime(1999,8,1,0,0,0);
time_sec = seconds(dates - t_ref);
Nt = length(time_sec);
nccreate(input_nc,'time', ...
    'Dimensions',{'time',Nt}, ...
    'Datatype','double');
ncwrite(input_nc,'time',time_sec);

%*** writing coordinate of wave nodes
Ns = length(ysave); %number of nodes
crs = projcrs(2193); %coordinates need to be in meters
[station_x, station_y] = projfwd(crs, xsave, ysave);
nccreate(input_nc,'station_x', ...
    'Dimensions',{'station',Ns}, ...
    'Datatype','double');
nccreate(input_nc,'station_y', ...
    'Dimensions',{'station',Ns}, ...
    'Datatype','double');
ncwrite(input_nc,'station_x',station_x);
ncwrite(input_nc,'station_y',station_y);

%************interpolating and writing in the netcdf file 

%****interpolation weights
% The Oceanum node positions do not change in time, and linear interpolation
% (with nearest-neighbour extrapolation) is linear in the node values. So the
% interpolation from the nodes to the contour points can be written as a
% weight matrix W (contour points x nodes), built once by interpolating unit
% vectors. Every hourly field is then interpolated with a single matrix
% product, instead of evaluating the scatteredInterpolant for every hour.
% Hours with missing (NaN) node values are interpolated as before, one by one,
% because 0*NaN would otherwise spread the NaN to all contour points.
Nn = numel(lon);
F_w = scatteredInterpolant(lon, lat, zeros(Nn,1), 'linear', 'nearest');
W = zeros(Ns,Nn);
for k = 1:Nn
    e_k = zeros(Nn,1);
    e_k(k) = 1;
    F_w.Values = e_k;
    W(:,k) = F_w(ysave,xsave);
end

%****hs along countour
hs_int = W*hs_idx; %nodes vs. time
for t = find(any(isnan(hs_idx),1))
    F_w.Values = hs_idx(:,t);
    hs_int(:,t) = F_w(ysave,xsave);
end

hs_int_m = hs_int'; %time vs. nodes
nccreate(input_nc,'point_hm0', ...
    'Dimensions',{'time',Nt,'station',Ns}, ...
    'Datatype','double');
ncwrite(input_nc,'point_hm0',hs_int_m);

%***dpm along countour
u = cosd(dpm_idx);
v = sind(dpm_idx);

u_int = W*u;
v_int = W*v;
for it = find(any(isnan(u),1) | any(isnan(v),1))
    F_w.Values = u(:,it);
    u_int(:,it) = F_w(ysave,xsave);
    F_w.Values = v(:,it);
    v_int(:,it) = F_w(ysave,xsave);
end
dpm_int = mod(atan2d(v_int,u_int),360);
dpm_int_m = dpm_int';

nccreate(input_nc,'point_wavdir', ...
    'Dimensions',{'time',Nt,'station',Ns}, ...
    'Datatype','double');
ncwrite(input_nc,'point_wavdir',dpm_int_m);

%***wave period along countour
tps_int = W*tps_idx;
for pt = find(any(isnan(tps_idx),1))
    F_w.Values = tps_idx(:,pt);
    tps_int(:,pt) = F_w(ysave,xsave);
end
tps_int_m = tps_int';
nccreate(input_nc,'point_tp', ...
    'Dimensions',{'time',Nt,'station',Ns}, ...
    'Datatype','double');
ncwrite(input_nc,'point_tp',tps_int_m);

ncwriteatt(input_nc,'time','units','seconds since 1999-08-01 00:00:00');

ncwriteatt(input_nc,'time', 'long_name','time');

ncwriteatt(input_nc,'point_hm0', ...
    'units', 'm');

ncwriteatt(input_nc,'point_tp', ...
    'units', 's');

ncwriteatt(input_nc,'point_wavdir', ...
    'units', 'degrees');

%saving variables
data.dates      = dates;
data.xsave      = xsave;
data.ysave      = ysave;
data.dpm_int_m  = dpm_int_m;
data.hs_int_m   = hs_int_m;
data.tps_int_m  = tps_int_m;
save('wave_data.mat','data','-v7.3')






