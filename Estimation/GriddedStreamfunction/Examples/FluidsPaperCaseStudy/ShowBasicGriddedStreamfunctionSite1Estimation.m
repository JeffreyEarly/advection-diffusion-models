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

fit = GriddedStreamfunction(trajectories, psiS=[2 2 3]);

tDays = t/86400;
f0 = 2 * 7.2921e-5 * sin(siteData.lat0*pi/180);
mx = fit.centerOfMassTrajectory.x(t);
my = fit.centerOfMassTrajectory.y(t);
uObserved = reshape(fit.observedXVelocity, numel(t), nDrifters);
vObserved = reshape(fit.observedYVelocity, numel(t), nDrifters);
uMesoscale = reshape(fit.uMesoscaleObserved, numel(t), nDrifters);
vMesoscale = reshape(fit.vMesoscaleObserved, numel(t), nDrifters);
uSubmesoscale = reshape(fit.uSubmesoscaleObserved, numel(t), nDrifters);
vSubmesoscale = reshape(fit.vSubmesoscaleObserved, numel(t), nDrifters);
qObserved = reshape(fit.centeredX, numel(t), nDrifters);
rObserved = reshape(fit.centeredY, numel(t), nDrifters);
qMesoscale = zeros(size(x));
rMesoscale = zeros(size(y));
for iDrifter = 1:nDrifters
    qMesoscale(:, iDrifter) = fit.centeredMesoscaleTrajectories(iDrifter).x(t);
    rMesoscale(:, iDrifter) = fit.centeredMesoscaleTrajectories(iDrifter).y(t);
end

uBackground = fit.uBackground(t);
vBackground = fit.vBackground(t);
uBackgroundGrid = repmat(uBackground, 1, nDrifters);
vBackgroundGrid = repmat(vBackground, 1, nDrifters);
sigmaN = fit.sigma_n(t, mx, my);
sigmaS = fit.sigma_s(t, mx, my);
sigma = hypot(sigmaN, sigmaS);
zeta = fit.zeta(t, mx, my);

% Use the first drifter as the default one-trajectory example.
iTrajectory = 1;
uReconstruction = uMesoscale + uBackgroundGrid + uSubmesoscale;
vReconstruction = vMesoscale + vBackgroundGrid + vSubmesoscale;
uDecompositionError = max(abs(uObserved(:, iTrajectory) - uReconstruction(:, iTrajectory)));
vDecompositionError = max(abs(vObserved(:, iTrajectory) - vReconstruction(:, iTrajectory)));
if uDecompositionError > 1e-10 || vDecompositionError > 1e-10
    error("GriddedStreamfunction:DecompositionMismatch", ...
        "The raw decomposition failed to reconstruct the observed velocity for drifter %d.", iTrajectory);
end

fprintf("Site 1 gridded streamfunction fit\n");
fprintf("  drifters: %d\n", nDrifters);
fprintf("  samples:  %d\n", numel(t));
fprintf("  mean sigma/f0: %.4f\n", mean(sigma/f0));
fprintf("  mean zeta/f0:  %.4f\n", mean(zeta/f0));
fprintf("  mean u_bg:     %.4f m/s\n", mean(uBackground));
fprintf("  mean v_bg:     %.4f m/s\n", mean(vBackground));
fprintf("  drifter %d max |u-u_{recon}|: %.3e m/s\n", iTrajectory, uDecompositionError);
fprintf("  drifter %d max |v-v_{recon}|: %.3e m/s\n", iTrajectory, vDecompositionError);

figure(Color="w");
tl = tiledlayout(3, 1, TileSpacing="compact", Padding="compact");

ax1 = nexttile;
plot(ax1, tDays, sigma/f0, LineWidth=2);
hold(ax1, "on");
plot(ax1, tDays, zeta/f0, LineWidth=2);
ylabel(ax1, "\sigma/f_0, \zeta/f_0");
xlim(ax1, [tDays(1), tDays(end)]);
legend(ax1, "\sigma/f_0", "\zeta/f_0", Location="best");
box(ax1, "on");
title(ax1, "Site 1 gridded streamfunction fit");

ax2 = nexttile;
plot(ax2, tDays, uBackground, LineWidth=2);
hold(ax2, "on");
plot(ax2, tDays, vBackground, LineWidth=2);
ylabel(ax2, "u_{bg}, v_{bg} (m/s)");
xlim(ax2, [tDays(1), tDays(end)]);
legend(ax2, "u_{bg}", "v_{bg}", Location="best");
box(ax2, "on");

ax3 = nexttile;
plot(ax3, qObserved/1000, rObserved/1000, Color=0.8*[1 1 1]);
hold(ax3, "on");
plot(ax3, qMesoscale/1000, rMesoscale/1000, LineWidth=1.5);
axis(ax3, "equal");
xlabel(ax3, "q (km)");
ylabel(ax3, "r (km)");
legend(ax3, "Observed", "Mesoscale fit", Location="best");
box(ax3, "on");

title(tl, "GriddedStreamfunction Site 1");

figure(Color="w");
tlDecomposition = tiledlayout(2, 1, TileSpacing="compact", Padding="compact");

axU = nexttile;
hold(axU, "on");
scatter(axU, tDays, uObserved(:, iTrajectory), 5^2, "k", "filled", DisplayName="observed");
plot(axU, tDays, uReconstruction(:, iTrajectory), "k", LineWidth=1.5, DisplayName="reconstruction");
plot(axU, tDays, uMesoscale(:, iTrajectory), LineWidth=1.5, Color=[0 0.4470 0.7410], DisplayName="mesoscale");
plot(axU, tDays, uBackgroundGrid(:, iTrajectory), LineWidth=1.5, Color=[0.8500 0.3250 0.0980], DisplayName="background");
plot(axU, tDays, uSubmesoscale(:, iTrajectory), LineWidth=1.5, Color=[0.4660 0.6740 0.1880], DisplayName="submesoscale");
ylabel(axU, "u (m/s)");
xlim(axU, [tDays(1), tDays(end)]);
legend(axU, Location="best");
box(axU, "on");

axV = nexttile;
hold(axV, "on");
scatter(axV, tDays, vObserved(:, iTrajectory), 5^2, "k", "filled", DisplayName="observed");
plot(axV, tDays, vReconstruction(:, iTrajectory), "k", LineWidth=1.5, DisplayName="reconstruction");
plot(axV, tDays, vMesoscale(:, iTrajectory), LineWidth=1.5, Color=[0 0.4470 0.7410], DisplayName="mesoscale");
plot(axV, tDays, vBackgroundGrid(:, iTrajectory), LineWidth=1.5, Color=[0.8500 0.3250 0.0980], DisplayName="background");
plot(axV, tDays, vSubmesoscale(:, iTrajectory), LineWidth=1.5, Color=[0.4660 0.6740 0.1880], DisplayName="submesoscale");
xlabel(axV, "time (days)");
ylabel(axV, "v (m/s)");
xlim(axV, [tDays(1), tDays(end)]);
legend(axV, Location="best");
box(axV, "on");

title(tlDecomposition, sprintf("Site 1 drifter %d raw velocity decomposition", iTrajectory));
