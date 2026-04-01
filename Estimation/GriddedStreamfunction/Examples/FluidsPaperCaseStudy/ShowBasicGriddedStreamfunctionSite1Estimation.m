% Show a simple Site 1 fit using the current gridded-streamfunction estimator.
scriptDir = fileparts(mfilename("fullpath"));
dataPath = fullfile(scriptDir, "..", "..", "..", "ExampleData", "LatMix2011", "smoothedGriddedRho1Drifters.mat");
if ~isfile(dataPath)
    error("GriddedStreamfunction:MissingExampleData", "Expected local example data at %s.", dataPath);
end

siteData = load(dataPath);
t = reshape(siteData.t, [], 1);
x = siteData.x(:, 1:(end - 1));
y = siteData.y(:, 1:(end - 1));
nDrifters = size(x, 2);

trajectories = TrajectorySpline.empty(0, 1);
for iDrifter = 1:nDrifters
    trajectories(end + 1, 1) = TrajectorySpline(t, x(:, iDrifter), y(:, iDrifter), S=3);
end

fit = GriddedStreamfunction.fitFromTrajectorySplines(trajectories);

tDays = t/86400;
f0 = 2 * 7.2921e-5 * sin(siteData.lat0*pi/180);
mx = fit.centerOfMassTrajectory.x(t);
my = fit.centerOfMassTrajectory.y(t);
tGrid = repmat(t, 1, nDrifters);
mxGrid = repmat(mx, 1, nDrifters);
myGrid = repmat(my, 1, nDrifters);

uMesoscale = fit.uMesoscale(tGrid, x, y);
vMesoscale = fit.vMesoscale(tGrid, x, y);
qObserved = x - mxGrid;
rObserved = y - myGrid;
qMesoscale = qObserved(1, :) + cumtrapz(t, uMesoscale);
rMesoscale = rObserved(1, :) + cumtrapz(t, vMesoscale);

uBackground = fit.uBackground(t);
vBackground = fit.vBackground(t);
sigmaN = fit.sigma_n(t, mx, my);
sigmaS = fit.sigma_s(t, mx, my);
sigma = hypot(sigmaN, sigmaS);
zeta = fit.zeta(t, mx, my);

fprintf("Site 1 gridded streamfunction fit\n");
fprintf("  drifters: %d\n", nDrifters);
fprintf("  samples:  %d\n", numel(t));
fprintf("  mean sigma/f0: %.4f\n", mean(sigma/f0));
fprintf("  mean zeta/f0:  %.4f\n", mean(zeta/f0));
fprintf("  mean u_bg:     %.4f m/s\n", mean(uBackground));
fprintf("  mean v_bg:     %.4f m/s\n", mean(vBackground));

figure("Color", "w");
tl = tiledlayout(3, 1, "TileSpacing", "compact", "Padding", "compact");

ax1 = nexttile;
plot(ax1, tDays, sigma/f0, "LineWidth", 2);
hold(ax1, "on");
plot(ax1, tDays, zeta/f0, "LineWidth", 2);
ylabel(ax1, "\sigma/f_0, \zeta/f_0");
xlim(ax1, [tDays(1), tDays(end)]);
legend(ax1, "\sigma/f_0", "\zeta/f_0", "Location", "best");
box(ax1, "on");
title(ax1, "Site 1 gridded streamfunction fit");

ax2 = nexttile;
plot(ax2, tDays, uBackground, "LineWidth", 2);
hold(ax2, "on");
plot(ax2, tDays, vBackground, "LineWidth", 2);
ylabel(ax2, "u_{bg}, v_{bg} (m/s)");
xlim(ax2, [tDays(1), tDays(end)]);
legend(ax2, "u_{bg}", "v_{bg}", "Location", "best");
box(ax2, "on");

ax3 = nexttile;
plot(ax3, qObserved/1000, rObserved/1000, "Color", 0.8*[1 1 1]);
hold(ax3, "on");
plot(ax3, qMesoscale/1000, rMesoscale/1000, "LineWidth", 1.5);
axis(ax3, "equal");
xlabel(ax3, "q (km)");
ylabel(ax3, "r (km)");
legend(ax3, "Observed", "Mesoscale fit", "Location", "best");
box(ax3, "on");

title(tl, "GriddedStreamfunction Site 1");
