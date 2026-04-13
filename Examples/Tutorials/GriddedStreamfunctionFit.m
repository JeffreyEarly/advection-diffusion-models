%% Tutorial Metadata
% Title: Gridded streamfunction fit
% Slug: gridded-streamfunction-fit
% Description: Fit a zero-vorticity gridded streamfunction to the Site 1 LatMix drifter cluster and inspect the fitted decomposition.
% NavOrder: 2

%% Load the Site 1 drifters
% `GriddedStreamfunction` fits one common background path, one centered
% mesoscale streamfunction, and one submesoscale residual per drifter so
% that
%
% $$ \dot{x}_k = u^{\mathrm{meso}} + u^{\mathrm{bg}} + u_k^{\mathrm{sm}}, \qquad \dot{y}_k = v^{\mathrm{meso}} + v^{\mathrm{bg}} + v_k^{\mathrm{sm}}. $$
%
% The checked-in Site 1 LatMix cluster translates together in the fixed
% frame while also shearing and spreading relative to its center of mass.
% This tutorial starts from those observed trajectories, then removes the
% common translation to fit a zero-vorticity mesoscale model.
scriptDir = fileparts(mfilename("fullpath"));
dataPath = fullfile(scriptDir, "..", "..", "Estimation", "ExampleData", "LatMix2011", "smoothedGriddedRho1Drifters.mat");
if ~isfile(dataPath)
    error("GriddedStreamfunction:MissingExampleData", "Expected local example data at %s.", dataPath);
end

siteData = load(dataPath);
t = reshape(siteData.t, [], 1);
x = siteData.x(:, 1:(end - 1));
y = siteData.y(:, 1:(end - 1));
nDrifters = size(x, 2);
tDays = t/86400;
f0 = 2 * 7.2921e-5 * sin(siteData.lat0*pi/180);

%% Plot the drifters in the fixed frame
% In the laboratory frame, the cluster translation and the relative motion
% are superimposed. The fit will later separate those two pieces.
figure(Color="w", Position=[100 100 430 360]); axFixed = axes; hold(axFixed, "on")
for iDrifter = 1:nDrifters
    plot(axFixed, x(:, iDrifter)/1000, y(:, iDrifter)/1000, LineWidth=1.2);
end
axis(axFixed, "equal"); xlabel(axFixed, "x (km)"); ylabel(axFixed, "y (km)")
title(axFixed, "Fixed-frame drifters"); box(axFixed, "on")
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("fixed-frame-drifters", Caption="In the fixed frame, the Site 1 drifters translate together while their relative spreading remains embedded in the common cluster motion."); end

%% Convert the drifters to trajectory splines
% Each drifter becomes one cubic `TrajectorySpline`, so the estimator
% works with both positions and spline-derived velocities on the same
% support times.
trajectoryCell = cell(nDrifters, 1);
for iDrifter = 1:nDrifters
    trajectoryCell{iDrifter} = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=3);
end
trajectories = vertcat(trajectoryCell{:});

%% Fit a zero-vorticity mesoscale model
% The Site 1 tutorial uses a quadratic-in-space, cubic-in-time mesoscale
% spline basis together with `mesoscaleConstraint="zeroVorticity"`, which
% keeps the fitted mesoscale flow harmonic in the centered frame.
fit = GriddedStreamfunction.fromTrajectories(trajectories, psiS=[2 2 3], fastS=3, mesoscaleConstraint="zeroVorticity");
decomposition = fit.decomposition;

%% Move to the center-of-mass frame
% The fitted center-of-mass trajectory defines the centered coordinates
%
% $$ q_k = x_k - m_x(t), \qquad r_k = y_k - m_y(t), $$
%
% which remove the common translation and leave the relative motion seen
% by the mesoscale spline.
figure(Color="w", Position=[100 100 430 360]); axCentered = axes; hold(axCentered, "on")
for iDrifter = 1:nDrifters
    trajectory = fit.observedTrajectories(iDrifter); ti = trajectory.t;
    [~, qi, ri] = fit.centeredCoordinates(ti, trajectory.x(ti), trajectory.y(ti));
    plot(axCentered, qi/1000, ri/1000, LineWidth=1.2);
end
axis(axCentered, "equal"); xlabel(axCentered, "q (km)"); ylabel(axCentered, "r (km)")
title(axCentered, "Center-of-mass-frame drifters"); box(axCentered, "on")
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("centered-frame-drifters", Caption="Removing the fitted center-of-mass translation reveals the relative motion that the mesoscale streamfunction and residual decomposition must explain."); end

%% Evaluate diagnostics on the center-of-mass path
% Evaluating the fit on the center-of-mass path gives a compact summary of
% the recovered strain, background drift, and the near-zero mesoscale
% vorticity implied by the chosen constraint.
xCom = fit.centerOfMassTrajectory.x(t);
yCom = fit.centerOfMassTrajectory.y(t);
sigmaN = fit.sigma_n(t, xCom, yCom);
sigmaS = fit.sigma_s(t, xCom, yCom);
sigma = hypot(sigmaN, sigmaS);
zeta = fit.zeta(t, xCom, yCom);
thetaDegrees = GriddedStreamfunction.visualPrincipalStrainAngle(sigmaN, sigmaS);
uBackground = fit.uBackground(t);
vBackground = fit.vBackground(t);

%% Plot the fitted diagnostics
figure(Color="w", Position=[100 100 700 560]); tlDiagnostics = tiledlayout(3, 1, TileSpacing="compact", Padding="compact");

axRate = nexttile; hold(axRate, "on")
plot(axRate, tDays, sigma/f0, LineWidth=1.5); plot(axRate, tDays, zeta/f0, LineWidth=1.5)
ylabel(axRate, "rate / f_0"); xlim(axRate, [tDays(1), tDays(end)])
legend(axRate, "\sigma", "\zeta", Location="best"); box(axRate, "on"); title(axRate, "Zero-vorticity diagnostics")

axTheta = nexttile;
plot(axTheta, tDays, thetaDegrees, LineWidth=1.5)
ylabel(axTheta, "\theta (deg)"); xlim(axTheta, [tDays(1), tDays(end)]); box(axTheta, "on")

axBackground = nexttile; hold(axBackground, "on")
plot(axBackground, tDays, uBackground, LineWidth=1.5); plot(axBackground, tDays, vBackground, LineWidth=1.5)
xlabel(axBackground, "time (days)"); ylabel(axBackground, "u_bg, v_bg (m/s)")
xlim(axBackground, [tDays(1), tDays(end)]); legend(axBackground, "u_bg", "v_bg", Location="best"); box(axBackground, "on")

title(tlDiagnostics, "Fitted diagnostics")
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("fitted-diagnostics", Caption="The zero-vorticity fit retains a time-varying strain field and background drift while keeping the mesoscale relative vorticity near zero along the fitted center-of-mass path."); end

%% Plot the fixed-frame decomposition
% The fixed-frame decomposition stores a common background path together
% with one mesoscale and one submesoscale trajectory for each drifter.
backgroundX = fit.backgroundTrajectory.x(t);
backgroundY = fit.backgroundTrajectory.y(t);

figure(Color="w", Position=[100 100 1080 340]); tlDecomposition = tiledlayout(1, 3, TileSpacing="none", Padding="compact");
axBackgroundPath = nexttile;
plot(axBackgroundPath, backgroundX/1000, backgroundY/1000, "k", LineWidth=1.5)
axis(axBackgroundPath, "equal"); xlabel(axBackgroundPath, "x (km)"); ylabel(axBackgroundPath, "y (km)")
title(axBackgroundPath, "Common background path"); box(axBackgroundPath, "on")

axMesoscale = nexttile; hold(axMesoscale, "on")
for iDrifter = 1:nDrifters
    trajectory = fit.observedTrajectories(iDrifter); ti = trajectory.t; mesoscale = decomposition.fixedFrame.mesoscale(iDrifter);
    plot(axMesoscale, mesoscale.x(ti)/1000, mesoscale.y(ti)/1000, LineWidth=1.2)
end
axis(axMesoscale, "equal"); xMeso = xlim(axMesoscale); yMeso = ylim(axMesoscale);
xlim(axMesoscale, xMeso + 0.05 * max(diff(xMeso), 1) * [-1 1]); ylim(axMesoscale, yMeso + 0.05 * max(diff(yMeso), 1) * [-1 1])
xlabel(axMesoscale, "x (km)"); title(axMesoscale, "Fixed-frame mesoscale"); axMesoscale.YTickLabel = []; box(axMesoscale, "on")

axSubmesoscale = nexttile; hold(axSubmesoscale, "on")
for iDrifter = 1:nDrifters
    trajectory = fit.observedTrajectories(iDrifter); ti = trajectory.t; submesoscale = decomposition.fixedFrame.submesoscale(iDrifter);
    plot(axSubmesoscale, submesoscale.x(ti)/1000, submesoscale.y(ti)/1000, LineWidth=1.2)
end
axis(axSubmesoscale, "equal"); xSubmeso = xlim(axSubmesoscale); ySubmeso = ylim(axSubmesoscale);
xlim(axSubmesoscale, xSubmeso + 0.05 * max(diff(xSubmeso), 1) * [-1 1]); ylim(axSubmesoscale, ySubmeso + 0.05 * max(diff(ySubmeso), 1) * [-1 1])
xlabel(axSubmesoscale, "x (km)"); title(axSubmesoscale, "Fixed-frame submesoscale"); axSubmesoscale.YTickLabel = []; box(axSubmesoscale, "on")

title(tlDecomposition, "Trajectory decomposition")
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("trajectory-decomposition", Caption="The fitted decomposition separates one common translating background path from the coherent mesoscale motion and the smaller drifter-to-drifter residual excursions."); end

%% Reconstruct the velocity of one drifter
% For an individual drifter, the spline-derived velocity is reconstructed
% directly from the fitted component velocities,
%
% $$ \mathbf{u}^{\mathrm{obs}}_k = \mathbf{u}^{\mathrm{bg}} + \mathbf{u}^{\mathrm{meso}}_k + \mathbf{u}^{\mathrm{sm}}_k. $$
iDrifter = 1;
trajectory = fit.observedTrajectories(iDrifter); ti = trajectory.t; tDaysDrifter = ti/86400;
background = decomposition.fixedFrame.background(iDrifter); mesoscale = decomposition.fixedFrame.mesoscale(iDrifter); submesoscale = decomposition.fixedFrame.submesoscale(iDrifter);

uObserved = trajectory.u(ti);
vObserved = trajectory.v(ti);
uBackgroundDrifter = background.u(ti);
vBackgroundDrifter = background.v(ti);
uMesoscale = mesoscale.u(ti);
vMesoscale = mesoscale.v(ti);
uSubmesoscale = submesoscale.u(ti);
vSubmesoscale = submesoscale.v(ti);
uReconstruction = uBackgroundDrifter + uMesoscale + uSubmesoscale;
vReconstruction = vBackgroundDrifter + vMesoscale + vSubmesoscale;

%% Print a compact fit summary
fprintf("Site 1 zero-vorticity fit\n");
fprintf("  drifters: %d\n", nDrifters);
fprintf("  max |zeta/f0| on COM path: %.3e\n", max(abs(zeta/f0)));
fprintf("  drifter %d max |u-u_recon|: %.3e m/s\n", iDrifter, max(abs(uObserved - uReconstruction)));
fprintf("  drifter %d max |v-v_recon|: %.3e m/s\n", iDrifter, max(abs(vObserved - vReconstruction)));

%% Plot the velocity decomposition for one drifter
figure(Color="w", Position=[100 100 760 420]); tlVelocity = tiledlayout(2, 1, TileSpacing="compact", Padding="compact");

axU = nexttile; hold(axU, "on")
scatter(axU, tDaysDrifter, uObserved, 5^2, "k", "filled", DisplayName="observed")
plot(axU, tDaysDrifter, uReconstruction, "k", LineWidth=1.5, DisplayName="reconstruction")
plot(axU, tDaysDrifter, uMesoscale, LineWidth=1.5, Color=[0 0.4470 0.7410], DisplayName="mesoscale")
plot(axU, tDaysDrifter, uBackgroundDrifter, LineWidth=1.5, Color=[0.8500 0.3250 0.0980], DisplayName="background")
plot(axU, tDaysDrifter, uSubmesoscale, LineWidth=1.5, Color=[0.4660 0.6740 0.1880], DisplayName="submesoscale")
ylabel(axU, "u (m/s)"); xlim(axU, [tDays(1), tDays(end)]); legend(axU, Location="best"); box(axU, "on")

axV = nexttile; hold(axV, "on")
scatter(axV, tDaysDrifter, vObserved, 5^2, "k", "filled", DisplayName="observed")
plot(axV, tDaysDrifter, vReconstruction, "k", LineWidth=1.5, DisplayName="reconstruction")
plot(axV, tDaysDrifter, vMesoscale, LineWidth=1.5, Color=[0 0.4470 0.7410], DisplayName="mesoscale")
plot(axV, tDaysDrifter, vBackgroundDrifter, LineWidth=1.5, Color=[0.8500 0.3250 0.0980], DisplayName="background")
plot(axV, tDaysDrifter, vSubmesoscale, LineWidth=1.5, Color=[0.4660 0.6740 0.1880], DisplayName="submesoscale")
xlabel(axV, "time (days)"); ylabel(axV, "v (m/s)")
xlim(axV, [tDays(1), tDays(end)]); legend(axV, Location="best"); box(axV, "on")

title(tlVelocity, sprintf("Velocity decomposition for drifter %d", iDrifter))
if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("velocity-decomposition", Caption="For an individual drifter, the fitted background, mesoscale, and submesoscale velocities add back to the observed spline-derived velocity at the sample times."); end
