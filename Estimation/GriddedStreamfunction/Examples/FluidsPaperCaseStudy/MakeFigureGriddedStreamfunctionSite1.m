% Make simple Site 1 figures directly from a gridded streamfunction fit.
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
f0 = 2 * 7.2921e-5 * sin(siteData.lat0*pi/180);

trajectories = TrajectorySpline.empty(0, 1);
for iDrifter = 1:nDrifters
    trajectories(end + 1, 1) = TrajectorySpline(t, x(:, iDrifter), y(:, iDrifter), S=3);
end

fit = GriddedStreamfunction(trajectories, psiS=[2 2 3]);

tDays = t/86400;
mx = fit.centerOfMassTrajectory.x(t);
my = fit.centerOfMassTrajectory.y(t);
sigmaN = fit.sigma_n(t, mx, my);
sigmaS = fit.sigma_s(t, mx, my);
sigma = hypot(sigmaN, sigmaS);
zeta = fit.zeta(t, mx, my);
thetaDegrees = principalThetaDegrees(unwrap(atan2(sigmaS, sigmaN))/2 * 180/pi);
uBackground = fit.uBackground(t);
vBackground = fit.vBackground(t);
backgroundX = fit.backgroundTrajectory.x(t);
backgroundY = fit.backgroundTrajectory.y(t);
xSubmesoscaleCell = cell(nDrifters, 1);
ySubmesoscaleCell = cell(nDrifters, 1);

figure(Color="w");
tl = tiledlayout(3, 1, TileSpacing="compact", Padding="compact");

ax1 = nexttile;
plot(ax1, tDays, sigma/f0, LineWidth=1.5);
hold(ax1, "on");
plot(ax1, tDays, zeta/f0, LineWidth=1.5);
ylabel(ax1, "sigma (f)");
xlim(ax1, [tDays(1), tDays(end)]);
legend(ax1, "\sigma", "\zeta", Location="best");
box(ax1, "on");
title(ax1, "Site 1 gridded streamfunction fit");

ax2 = nexttile;
plot(ax2, tDays, thetaDegrees, LineWidth=1.5);
ylabel(ax2, "\theta (deg)");
xlim(ax2, [tDays(1), tDays(end)]);
box(ax2, "on");

ax3 = nexttile;
plot(ax3, tDays, uBackground, LineWidth=1.5);
hold(ax3, "on");
plot(ax3, tDays, vBackground, LineWidth=1.5);
xlabel(ax3, "time (days)");
ylabel(ax3, "u_{bg}, v_{bg} (m/s)");
xlim(ax3, [tDays(1), tDays(end)]);
legend(ax3, "u_{bg}", "v_{bg}", Location="best");
box(ax3, "on");

title(tl, "Site 1 parameter time series");

figure(Color="w");
axMeso = axes;
hold(axMeso, "on");
for iDrifter = 1:nDrifters
    trajectory = fit.observedTrajectories(iDrifter);
    ti = trajectory.t;
    xMesoscale = fit.decomposition.fixedFrame.mesoscale(iDrifter).x(ti);
    yMesoscale = fit.decomposition.fixedFrame.mesoscale(iDrifter).y(ti);
    plot(axMeso, xMesoscale/1000, yMesoscale/1000, LineWidth=1.2);
end
axis(axMeso, "equal");
xlabel(axMeso, "x (km)");
ylabel(axMeso, "y (km)");
title(axMeso, "Site 1 fixed-frame mesoscale trajectories");
box(axMeso, "on");

figure(Color="w");
tlDecomp = tiledlayout(2, 1, TileSpacing="compact", Padding="compact");

axBackground = nexttile;
plot(axBackground, backgroundX/1000, backgroundY/1000, "k", LineWidth=1.5);
axis(axBackground, "equal");
ylabel(axBackground, "y (km)");
title(axBackground, "Common background path");
box(axBackground, "on");

axSubmeso = nexttile;
hold(axSubmeso, "on");
for iDrifter = 1:nDrifters
    trajectory = fit.observedTrajectories(iDrifter);
    ti = trajectory.t;
    xSubmesoscale = fit.decomposition.fixedFrame.submesoscale(iDrifter).x(ti);
    ySubmesoscale = fit.decomposition.fixedFrame.submesoscale(iDrifter).y(ti);
    xSubmesoscaleCell{iDrifter} = xSubmesoscale;
    ySubmesoscaleCell{iDrifter} = ySubmesoscale;
    plot(axSubmeso, xSubmesoscale/1000, ySubmesoscale/1000, LineWidth=1.2);
end
axis(axSubmeso, "equal");
xlabel(axSubmeso, "x (km)");
ylabel(axSubmeso, "y (km)");
title(axSubmeso, "Fixed-frame submesoscale trajectories");
box(axSubmeso, "on");

[decompositionXLimits, decompositionYLimits] = sharedTrajectoryLimits( ...
    [{backgroundX}; xSubmesoscaleCell], [{backgroundY}; ySubmesoscaleCell], 1000);
xlim(axBackground, decompositionXLimits)
xlim(axSubmeso, decompositionXLimits)
ylim(axBackground, decompositionYLimits)
ylim(axSubmeso, decompositionYLimits)

title(tlDecomp, "Site 1 fixed-frame decomposition");

figure(Color="w");
axCom = axes;
hold(axCom, "on");
for iDrifter = 1:nDrifters
    trajectory = fit.observedTrajectories(iDrifter);
    ti = trajectory.t;
    xi = trajectory.x(ti);
    yi = trajectory.y(ti);
    [~, qObserved, rObserved] = fit.centeredCoordinates(ti, xi, yi);
    qMesoscale = fit.decomposition.centeredFrame.mesoscale(iDrifter).x(ti);
    rMesoscale = fit.decomposition.centeredFrame.mesoscale(iDrifter).y(ti);

    if iDrifter == 1
        plot(axCom, qObserved/1000, rObserved/1000, Color=0.8*[1 1 1], DisplayName="Observed");
        plot(axCom, qMesoscale/1000, rMesoscale/1000, LineWidth=1.2, DisplayName="Mesoscale fit");
    else
        plot(axCom, qObserved/1000, rObserved/1000, Color=0.8*[1 1 1]);
        plot(axCom, qMesoscale/1000, rMesoscale/1000, LineWidth=1.2);
    end
end
axis(axCom, "equal");
xlabel(axCom, "q (km)");
ylabel(axCom, "r (km)");
title(axCom, "Site 1 COM-frame mesoscale trajectories");
legend(axCom, "Observed", "Mesoscale fit", Location="best");
box(axCom, "on");

function [xLimits, yLimits] = sharedTrajectoryLimits(xCell, yCell, scale)
xValues = cellfun(@(values) values(:), xCell, UniformOutput=false);
yValues = cellfun(@(values) values(:), yCell, UniformOutput=false);
xData = vertcat(xValues{:})/scale;
yData = vertcat(yValues{:})/scale;
xLimits = paddedLimits(xData);
yLimits = paddedLimits(yData);
end

function limits = paddedLimits(values)
limits = [min(values), max(values)];
range = diff(limits);
if range <= 0
    range = max(abs(limits(1)), 1);
end
limits = limits + 0.05 * range * [-1 1];
end

function thetaDegrees = principalThetaDegrees(thetaDegrees)
thetaDegrees = thetaDegrees(:);
thetaDegrees = mod(thetaDegrees + 45, 90) - 45;
end
