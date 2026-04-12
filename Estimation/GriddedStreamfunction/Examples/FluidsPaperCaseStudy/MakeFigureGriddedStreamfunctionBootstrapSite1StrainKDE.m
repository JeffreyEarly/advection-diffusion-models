scriptDir = fileparts(mfilename("fullpath"));
siteNumber = 1;
siteLabel = "Site 1";
dataFilename = "smoothedGriddedRho1Drifters.mat";
nBootstraps = 500;
randomSeed = 0;
scoreStride = 6;
psiS = [2 2 1];
fastS = 3;
mesoscaleConstraint = "zeroVorticity";
timeDays = [3 4 5];
shouldSaveFigures = 0;

figureFont = "Helvetica";
figureWidth = 700;
figureHeight = 220 + 180 * ceil(numel(timeDays)/2);
titleSize = 11;
axisLabelSize = 10;
tickSize = 9;
annotationSize = 9;
pointColor = [0 0.4470 0.7410];
ringColor = 0.75 * [1 1 1];
pctTarget = (0.8:-0.1:0.1).';
minimumHalfWidth = 1e-6;

dataPath = fullfile(scriptDir, "..", "..", "..", "ExampleData", "LatMix2011", dataFilename);
cacheDirectory = fullfile(scriptDir, "BootstrapData");
if exist(cacheDirectory, "dir") == 0
    mkdir(cacheDirectory);
end

psiSTag = join(string(psiS), "-");
cacheFilename = "Rho" + siteNumber + ...
    "GriddedStreamfunctionBootstrapFits" + nBootstraps + ...
    "_seed" + randomSeed + ...
    "_stride" + scoreStride + ...
    "_fastS" + fastS + ...
    "_psiS" + psiSTag + ...
    "_mesoscaleConstraint-" + mesoscaleConstraint + ".nc";
cachePath = fullfile(cacheDirectory, cacheFilename);

if isfile(cachePath)
    bootstrap = GriddedStreamfunctionBootstrap.fromFile(cachePath);
else
    siteData = load(dataPath);
    t = reshape(siteData.t, [], 1);
    x = siteData.x(:, 1:(end - 1));
    y = siteData.y(:, 1:(end - 1));

    trajectories = TrajectorySpline.empty(0, 1);
    for iDrifter = 1:size(x, 2)
        trajectories(end + 1, 1) = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=3); %#ok<SAGROW>
    end

    bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
        trajectories, ...
        nBootstraps=nBootstraps, ...
        randomSeed=randomSeed, ...
        queryTimes=t, ...
        scoreStride=scoreStride, ...
        psiS=psiS, ...
        fastS=fastS, ...
        mesoscaleConstraint=mesoscaleConstraint);
    bootstrap.writeToFile(cachePath, shouldOverwriteExisting=true);
end

timeDays = reshape(timeDays, [], 1);
tBootstrapDays = bootstrap.queryTimes/86400;
if any(timeDays < tBootstrapDays(1) | timeDays > tBootstrapDays(end))
    error("GriddedStreamfunction:TimeOutOfRange", ...
        "Requested timeDays must lie inside [%.6g, %.6g] days.", tBootstrapDays(1), tBootstrapDays(end));
end

timeIndices = zeros(size(timeDays));
for iTime = 1:numel(timeDays)
    [~, timeIndices(iTime)] = min(abs(tBootstrapDays - timeDays(iTime)));
end

if numel(unique(timeIndices)) ~= numel(timeIndices)
    error("GriddedStreamfunction:DuplicateTimeIndex", ...
        "Requested timeDays must map to distinct query samples.");
end

matchedTimeDays = tBootstrapDays(timeIndices);
disp(table(timeDays, matchedTimeDays, timeIndices, VariableNames=["RequestedDay", "MatchedDay", "QueryIndex"]))

sigmaNSelection = bootstrap.summary.sigma_n(timeIndices, :);
sigmaSSelection = bootstrap.summary.sigma_s(timeIndices, :);
fullSigmaN = bootstrap.fullSummary.sigma_n(timeIndices);
fullSigmaS = bootstrap.fullSummary.sigma_s(timeIndices);
allCoordinates = [sigmaNSelection(:); sigmaSSelection(:); fullSigmaN(:); fullSigmaS(:)];
halfWidth = 1.1 * max(abs(allCoordinates), [], "all");
if ~isfinite(halfWidth) || halfWidth <= 0
    halfWidth = minimumHalfWidth;
end
halfWidth = max(halfWidth, minimumHalfWidth);
minimum = [-halfWidth, -halfWidth];
maximum = [halfWidth, halfWidth];

figure(Color="w", Units="points", Position=[50 50 figureWidth figureHeight]);
tl = tiledlayout("flow", TileSpacing="compact", Padding="compact");
cmap = parula(256);
cmap(1, :) = 1;
colormap(cmap);

for iTime = 1:numel(timeDays)
    ax = nexttile(tl);
    samples = [ ...
        reshape(bootstrap.summary.sigma_n(timeIndices(iTime), :), [], 1), ...
        reshape(bootstrap.summary.sigma_s(timeIndices(iTime), :), [], 1)];
    fullPoint = [fullSigmaN(iTime), fullSigmaS(iTime)];

    densityModel = KernelDensityEstimate.fromData(samples, minimum=minimum, maximum=maximum);
    [density, gridVectors] = densityModel.densityOnGrid();
    contourLevels = DensityLevelForCDF(gridVectors, density, pctTarget);
    contourLevels = reshape(contourLevels, [], 1);
    contourLevels = contourLevels(isfinite(contourLevels));
    contourLevels = unique(contourLevels, "sorted");
    if numel(contourLevels) < 2
        densityMinimum = min(density, [], "all");
        densityMaximum = max(density, [], "all");
        if densityMaximum <= densityMinimum
            contourLevels = [densityMinimum; densityMinimum + eps(densityMinimum + 1)];
        else
            contourLevels = linspace(densityMinimum, densityMaximum, max(2, numel(pctTarget))).';
        end
    end

    hold(ax, "on");
    for radius = (2:2:10) * 1e-6
        rectangle(ax, Position=[-radius -radius 2 * radius 2 * radius], Curvature=[1 1], ...
            EdgeColor=ringColor, LineStyle=":", LineWidth=0.8);
    end
    contourf(ax, gridVectors{1}, gridVectors{2}, density.', contourLevels, LineStyle="none");
    scatter(ax, samples(:, 1), samples(:, 2), 12, pointColor, "filled", ...
        MarkerEdgeColor="none", MarkerFaceAlpha=0.18);
    scatter(ax, fullPoint(1), fullPoint(2), 42, "w", "filled", MarkerEdgeColor="k", LineWidth=1.1);
    axis(ax, "equal");
    xlim(ax, minimum(1) + [0, 1] * (maximum(1) - minimum(1)));
    ylim(ax, minimum(2) + [0, 1] * (maximum(2) - minimum(2)));
    box(ax, "on");
    xlabel(ax, "\sigma_n (s^{-1})", FontSize=axisLabelSize, FontName=figureFont);
    ylabel(ax, "\sigma_s (s^{-1})", FontSize=axisLabelSize, FontName=figureFont);
    title(ax, sprintf("t = %.3f d", matchedTimeDays(iTime)), FontSize=titleSize, FontName=figureFont);
    text(ax, 0.03, 0.97, sprintf("requested %.3f d\nindex %d", timeDays(iTime), timeIndices(iTime)), ...
        Units="normalized", HorizontalAlignment="left", VerticalAlignment="top", ...
        FontSize=annotationSize, FontName=figureFont, Color=[0.2 0.2 0.2]);
    set(ax, FontSize=tickSize, FontName=figureFont);
end

title(tl, siteLabel + " bootstrap strain KDE", FontSize=titleSize, FontName=figureFont);

if shouldSaveFigures == 1
    print(fullfile(scriptDir, "Site1BootstrapStrainKDE.eps"), "-depsc2");
end
