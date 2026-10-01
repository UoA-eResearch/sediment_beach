

function [s, t1, lat_s, lon_s, tit, sea_long, sea_lat] = load_data2(idx,site)

%This function load and  reprojects the coastsat coordinates to UTM and filters the
%years to input as initial condition for the sediment transport model.

%Creates a table for the initial/final data (first observation of 1999/last 2024 due to hindcast wave data) and saves
%the output files (just initial data) as txt files according to shoreline s input seetings

%after a correction is applied to ensure at least 3 points whithin
%a coastal section (separated by NaN)

if site == 1
    
    %files from north to south
    files  = {'nzd0299','nzd0300','nzd0301','nzd0303','nzd0307', ...
    'nzd0309','nzd0311','nzd0315','nzd0321','nzd0325','nzd0332'};
    fname = files{idx} + ".xlsx";
    fname=append(fname);
elseif site == 2
  
%files from west to east
    files  = {'nzd0204','nzd0207','nzd0217','nzd0222','nzd0220','nzd0226','nzd0229','nzd0231','nzd0234','nzd0239','nzd0240','nzd0238', 'nzd0236', 'nzd0233', 'nzd0230', 'nzd0227'};
    fname = files{idx} + ".xlsx";
    fname=append(fname);  
end

for i = 1:length(fname)
    data = readtable(input_file, 'Sheet', 'Intersect points'); %read coordinates of shoreline position
    data.dates = datetime(data.dates, 'InputFormat', 'yyyy-MM-dd HH:mm:ssXXX', 'TimeZone', 'UTC'); %convert dates
    utm_proj = projcrs(32760);
    data_utm = data(:, 1)
end