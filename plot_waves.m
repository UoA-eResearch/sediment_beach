function plot_waves

a = load('wave_data.mat');
t = a.data.dates;
ll= a.data.ysave;
h = a.data.hs_int_m;
w = a.data.dpm_int_m;
p = a.data.tps_int_m;

figure
hold on
% ---------------- Hs ----------------
subplot(3,1,1)
pcolor(t, ll, w');shading flat
set(gca,'YDir','reverse')
ylabel('Longitude')
ax = gca;
ax.XTickLabel = [];
colormap(flipud(cmocean('phase',5)))
cb2 = colorbar;
cb2 = colorbar;
cb2.Ticks = [0.1 90 180 270 350];
cb2.TickLabels = {'0','90','180','270','360'};
ylabel(cb2,'Wave direction (°)')



% ---------------- Wave direction ----------------
subplot(3,1,2)

pcolor(t, ll, h');shading flat
set(gca,'YDir','reverse')

ylabel('Longitude')
ax = gca;
ax.XTickLabel = [];
cmocean('haline');
cb1 = colorbar;
ylabel(cb1,'H_s (m)')


% ---------------- Wave period ----------------
subplot(3,1,3)

pcolor(t, ll, p')
set(gca,'YDir','reverse')
ylabel('Longitude')
ax = gca;
ax.XAxis.TickLabelFormat = 'yyyy';
ax.XTickLabel = string(year(ax.XTick));
shading flat
cmocean('thermal')
cb3 = colorbar;
ylabel(cb3,'T_p (s)')
xlabel('Timestep = 3 hours')
hold off




