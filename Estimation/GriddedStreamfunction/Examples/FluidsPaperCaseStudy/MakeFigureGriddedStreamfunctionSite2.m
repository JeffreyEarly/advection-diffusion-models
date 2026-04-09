scriptDir = fileparts(mfilename('fullpath'));
dataPath = fullfile(scriptDir, '..', '..', '..', 'ExampleData', 'LatMix2011', 'smoothedGriddedRho2Drifters.mat');
if ~isfile(dataPath)
    error("GriddedStreamfunction:MissingExampleData", "Expected local example data at %s.", dataPath);
end

siteData = load(dataPath);
t = reshape(siteData.t, [], 1);
x = siteData.x(:, 1:(end - 1));
y = siteData.y(:, 1:(end - 1));
nDrifters = size(x, 2);
[timeBasisMatrix, tKnotPoints, timeSplineDegree] = legacyTimeBasis(t, 4);
psiS = [2 2 timeSplineDegree];
mxLegacy = timeBasisMatrix * (timeBasisMatrix \ mean(x, 2));
myLegacy = timeBasisMatrix * (timeBasisMatrix \ mean(y, 2));
psiKnotPoints = { ...
    paddedQuadraticDomain(x - mxLegacy) ...
    paddedQuadraticDomain(y - myLegacy) ...
    tKnotPoints ...
    };

trajectories = TrajectorySpline.empty(0, 1);
for iDrifter = 1:nDrifters
    xSpline = ConstrainedSpline(t, x(:, iDrifter), S=3, knotPoints=tKnotPoints);
    ySpline = ConstrainedSpline(t, y(:, iDrifter), S=3, knotPoints=tKnotPoints);
    trajectories(end + 1, 1) = TrajectorySpline.fromComponentSplines(t, xSpline, ySpline); %#ok<AGROW>
end

fit = GriddedStreamfunction(trajectories, psiKnotPoints=psiKnotPoints, ...
    psiS=psiS, fastKnotPoints=tKnotPoints, fastS=timeSplineDegree);

scaleFactor = 1;
run(fullfile(scriptDir, 'LoadFigureDefaults.m'));

tDays = t/86400;
f0 = 2 * 7.2921e-5 * sin(siteData.lat0*pi/180);
mx = fit.centerOfMassTrajectory.x(t);
my = fit.centerOfMassTrajectory.y(t);
sigmaN = fit.sigma_n(t, mx, my);
sigmaS = fit.sigma_s(t, mx, my);
sigma = hypot(sigmaN, sigmaS);
zeta = fit.zeta(t, mx, my);
thetaDegrees = unwrap(atan2(sigmaS, sigmaN))/2 * 180/pi;
uBackground = fit.uBackground(t);
vBackground = fit.vBackground(t);
uObserved = zeros(numel(t), nDrifters);
vObserved = zeros(numel(t), nDrifters);
uMesoscale = zeros(numel(t), nDrifters);
vMesoscale = zeros(numel(t), nDrifters);
qObservedCell = cell(nDrifters, 1);
rObservedCell = cell(nDrifters, 1);
qMesoscaleCell = cell(nDrifters, 1);
rMesoscaleCell = cell(nDrifters, 1);
for iDrifter = 1:nDrifters
    trajectory = fit.observedTrajectories(iDrifter);
    ti = trajectory.t;
    xi = trajectory.x(ti);
    yi = trajectory.y(ti);
    [~, qObserved, rObserved] = fit.centeredCoordinates(ti, xi, yi);
    qMesoscale = fit.decomposition.centeredFrame.mesoscale(iDrifter).x(ti);
    rMesoscale = fit.decomposition.centeredFrame.mesoscale(iDrifter).y(ti);

    uObserved(:, iDrifter) = trajectory.u(ti);
    vObserved(:, iDrifter) = trajectory.v(ti);
    uMesoscale(:, iDrifter) = fit.uMesoscale(ti, xi, yi);
    vMesoscale(:, iDrifter) = fit.vMesoscale(ti, xi, yi);
    qObservedCell{iDrifter} = qObserved;
    rObservedCell{iDrifter} = rObserved;
    qMesoscaleCell{iDrifter} = qMesoscale;
    rMesoscaleCell{iDrifter} = rMesoscale;
end
uBackgroundGrid = repmat(uBackground, 1, nDrifters);
vBackgroundGrid = repmat(vBackground, 1, nDrifters);
uSubmesoscale = uObserved - uMesoscale - uBackgroundGrid;
vSubmesoscale = vObserved - vMesoscale - vBackgroundGrid;
observedSpeedRms = sqrt(mean(uObserved.^2 + vObserved.^2, 2));
mesoscaleSpeedRms = sqrt(mean(uMesoscale.^2 + vMesoscale.^2, 2));
backgroundSpeedRms = sqrt(uBackground.^2 + vBackground.^2);
submesoscaleSpeedRms = sqrt(mean(uSubmesoscale.^2 + vSubmesoscale.^2, 2));

figure('Units', 'points', 'Position', [50 50 figureWidthMedium 520*scaleFactor])
set(gcf, 'PaperPositionMode', 'auto')
set(gcf, 'Color', 'w');

tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
plot(tDays, zeta/f0, 'LineWidth', 2)
ylabel('\zeta / f_0', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
xlim([tDays(1) tDays(end)])
set(gca, 'FontSize', figureAxisTickSize, 'FontName', figureFont)
title('Site 2 Gridded Streamfunction Fit', 'FontSize', figureTitleSize, 'FontName', figureFont)

nexttile
plot(tDays, sigma/f0, 'LineWidth', 2)
ylabel('\sigma / f_0', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
xlim([tDays(1) tDays(end)])
set(gca, 'FontSize', figureAxisTickSize, 'FontName', figureFont)

nexttile
plot(tDays, thetaDegrees, 'LineWidth', 2)
ylabel('\theta (°)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
xlim([tDays(1) tDays(end)])
set(gca, 'FontSize', figureAxisTickSize, 'FontName', figureFont)

nexttile
plot(tDays, uBackground, 'LineWidth', 2)
hold on
plot(tDays, vBackground, 'LineWidth', 2)
xlabel('time (days)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
ylabel('u_{bg}, v_{bg} (m/s)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
xlim([tDays(1) tDays(end)])
legend({'u_{bg}', 'v_{bg}'}, 'Location', 'best', 'NumColumns', 2, 'FontSize', figureLegendSize)
set(gca, 'FontSize', figureAxisTickSize, 'FontName', figureFont)

figure('Units', 'points', 'Position', [50 50 figureWidthOneColumn 420*scaleFactor])
set(gcf, 'PaperPositionMode', 'auto')
set(gcf, 'Color', 'w');

tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
hold on
for iDrifter = 1:nDrifters
    if iDrifter == 1
        plot(qObservedCell{iDrifter}/1000, rObservedCell{iDrifter}/1000, 'Color', 0.8*[1 1 1], 'DisplayName', 'Observed')
        plot(qMesoscaleCell{iDrifter}/1000, rMesoscaleCell{iDrifter}/1000, 'LineWidth', 1.5, 'DisplayName', 'Mesoscale fit')
    else
        plot(qObservedCell{iDrifter}/1000, rObservedCell{iDrifter}/1000, 'Color', 0.8*[1 1 1])
        plot(qMesoscaleCell{iDrifter}/1000, rMesoscaleCell{iDrifter}/1000, 'LineWidth', 1.5)
    end
end
axis equal
xlabel('q (km)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
ylabel('r (km)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
title('Centered trajectories', 'FontSize', figureTitleSize, 'FontName', figureFont)
legend({'Observed', 'Mesoscale fit'}, 'Location', 'best', 'FontSize', figureLegendSize)
set(gca, 'FontSize', figureAxisTickSize, 'FontName', figureFont)

nexttile
plot(tDays, observedSpeedRms, 'Color', 0.2*[1 1 1], 'LineWidth', 1.5)
hold on
plot(tDays, mesoscaleSpeedRms, 'LineWidth', 2)
plot(tDays, backgroundSpeedRms, 'LineWidth', 2)
plot(tDays, submesoscaleSpeedRms, 'LineWidth', 2)
xlabel('time (days)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
ylabel('RMS speed (m/s)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
xlim([tDays(1) tDays(end)])
legend({'Observed', 'Mesoscale', 'Background', 'Submesoscale'}, 'Location', 'best', 'FontSize', figureLegendSize)
set(gca, 'FontSize', figureAxisTickSize, 'FontName', figureFont)

function [basisMatrix, knotPoints, splineDegree] = legacyTimeBasis(t, dof)
if dof == 1
    splineDegree = 0;
    knotPoints = [t(1); t(end)];
    basisMatrix = ones(numel(t), 1);
    return
end

splineDegree = min([dof 4]) - 1;
tData = linspace(t(1), t(end), dof + 1).';
timeData = (tData(1:(end - 1)) + tData(2:end))/2;
knotPoints = BSpline.knotPointsForDataPoints(timeData, S=splineDegree);
knotPoints(1:(splineDegree + 1)) = tData(1);
knotPoints((end - splineDegree):end) = tData(end);
basisMatrix = BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=splineDegree);
end

function knotPoints = paddedQuadraticDomain(values)
extent = max(abs(values), [], "all");
padding = max(50, 0.25 * max(extent, eps));
lowerBound = -(extent + padding);
upperBound = extent + padding;
knotPoints = [repmat(lowerBound, 3, 1); repmat(upperBound, 3, 1)];
end
