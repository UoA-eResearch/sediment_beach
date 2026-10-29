output_file = fullfile('C:\Users\mpul348\Documents\sediment_beach\python\output\28_oct_per_per', 'output.mat');
if ~exist(output_file, 'file')
    output_file = fullfile(pwd, file_out, 'output.mat');
end
raw_out = load(output_file);


PHI = raw_out.O.PHI; %nearshore waves angle in relation to north, clockwise increasing
dPHI = raw_out.O.dPHI;%nearshore waves angle in relation to the coastline
QS = raw_out.O.QS; %transport
dPHIbr = raw_out.O.dPHIbr;%wave angle at the point of breaking in relation to the coastline
PHIc = raw_out.O.PHIc; %coastline orientation calculated from the coastline geometry.
PHIf = raw_out.O.PHIf;%smoothed foreshore orientation used to calculate the wave angle relative to the coast.
%PHIcs: coastline orientation used for the refraction calculation in the code you shared.


figure
tiledlayout(3,1,'TileSpacing','compact','Padding','compact')

nexttile %nearshore waves
% mean_PHI = atan2d(mean(sind(PHI), 'all', 'omitnan'),mean(cosd(PHI), 'all', 'omitnan'));
% PHI_anom = mod(PHI - mean_PHI + 180, 360) - 180;
pcolor(PHI_anom)
shading flat
colorbar
colormap(cmocean('phase'))
clim([-90 90])

nexttile %nearshore waves in relation to the coastline
pcolor(dPHI)
shading flat
colorbar
colormap(cmocean('balance'))
clim([-90 90])


nexttile %transport
pcolor(QS)
shading flat
colorbar
colormap(cmocean('balance'))
clim([-max(abs(QS(:)), [], 'omitnan'), ...
       max(abs(QS(:)), [], 'omitnan')])



disp(table( ...
    PHI(1:20)', ...
    PHIf(1:20)', ...
    dPHIbr(1:20)', ...
    'VariableNames', {'PHI','PHIf','dPHIbr'}))

figure;pcolor(PHIf);shading flat