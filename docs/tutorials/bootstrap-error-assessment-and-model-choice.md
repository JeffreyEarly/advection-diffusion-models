---
layout: default
title: Bootstrap error assessment and model choice
parent: Tutorials
nav_order: 3
mathjax: true
permalink: /tutorials/bootstrap-error-assessment-and-model-choice
---

# Bootstrap error assessment and model choice

Start from the deterministic gridded streamfunction fit, then bootstrap Site 1 to assess uncertainty and compare a few nearby model choices.

Source: `Examples/Tutorials/BootstrapErrorAssessmentAndModelChoice.m`

## Bootstrap the Site 1 fit

Continue from [Gridded streamfunction fit](./gridded-streamfunction-fit):
that tutorial builds the deterministic Site 1 zero-vorticity fit, and
this follow-on tutorial bootstraps the same workflow to assess
uncertainty and compare a few nearby model choices.

```matlab
scriptDir = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(scriptDir, "..", "..", "..", "jLab")));
siteData = load(fullfile(scriptDir, "..", "..", "Estimation", "ExampleData", "LatMix2011", "smoothedGriddedRho1Drifters.mat"));
t = reshape(siteData.t, [], 1);
x = siteData.x(:, 1:(end - 1));
y = siteData.y(:, 1:(end - 1));
tDays = t/86400;
f0 = 2 * 7.2921e-5 * sin(siteData.lat0*pi/180);

trajectoryCell = cell(size(x, 2), 1);
for iDrifter = 1:size(x, 2)
    trajectoryCell{iDrifter} = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=3);
end
trajectories = vertcat(trajectoryCell{:});

bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
    trajectories, nBootstraps=200, randomSeed=0, queryTimes=t, ...
    scoreStride=6, psiS=[2 2 1], fastS=3, mesoscaleConstraint="zeroVorticity");
restartPath = tempname + ".nc";
cleanupRestart = onCleanup(@() delete(restartPath)); %#ok<NASGU>
restartFile = bootstrap.writeToFile(restartPath, shouldOverwriteExisting=true);
restartFile.close();
```

## View one bootstrap strain KDE

A single day-4 slice shows what the whole-drifter bootstrap has done:
each point is one bootstrap fit in `(\sigma_n,\sigma_s)` space, and the
black marker is the full-data fit.

```matlab
kdeDay = 4;
[~, iKDE] = min(abs(tDays - kdeDay));
sigmaN = reshape(bootstrap.summary.sigma_n(iKDE, :), [], 1)/f0;
sigmaS = reshape(bootstrap.summary.sigma_s(iKDE, :), [], 1)/f0;
fullSigmaN = bootstrap.fullSummary.sigma_n(iKDE)/f0;
fullSigmaS = bootstrap.fullSummary.sigma_s(iKDE)/f0;
halfWidth = max(0.06, 1.1 * max(abs([sigmaN; sigmaS; fullSigmaN; fullSigmaS])));
ringRadii = (2:2:10) * 1e-6 / f0;
kdeStatistics = KernelDensityEstimate.planarStatisticsFromData( ...
    [sigmaN, sigmaS], referencePoint=[fullSigmaN, fullSigmaS], ...
    minimum=[-halfWidth, -halfWidth], maximum=[halfWidth, halfWidth], summaryMass=0.68);

figure(Color="w", Position=[100 100 430 360]);
axSlice = axes;
cmap = parula(256);
cmap(1, :) = 1;
colormap(cmap);
KernelDensityEstimate.plotPlanarStatistics( ...
    axSlice, kdeStatistics, samples=[sigmaN, sigmaS], contourMasses=(0.8:-0.1:0.1).', ...
    showWedge=false, scatterColor=[0 0.4470 0.7410], ringRadii=ringRadii, ...
    ringColor=0.75 * [1 1 1], ringLineStyle=":", ringLineWidth=0.8);
hold(axSlice, "on")
plot(axSlice, fullSigmaN, fullSigmaS, "kp", MarkerSize=10, MarkerFaceColor="k")
axis(axSlice, "equal")
xlim(axSlice, halfWidth * [-1 1])
ylim(axSlice, halfWidth * [-1 1])
xlabel(axSlice, "\sigma_n / f_0")
ylabel(axSlice, "\sigma_s / f_0")
title(axSlice, sprintf("day %.1f bootstrap KDE", tDays(iKDE)))
box(axSlice, "on")
```

![At day 4, the whole-drifter bootstrap cloud in `(\sigma_n,\sigma_s)` space shows the spread of zero-vorticity fits around the full-data solution.](./bootstrap-error-assessment-and-model-choice/bootstrap-kde-slice.png)

*At day 4, the whole-drifter bootstrap cloud in `(\sigma_n,\sigma_s)` space shows the spread of zero-vorticity fits around the full-data solution.*

## Plot the bootstrap strain KDEs

A few matched KDE slices show how the cloud moves through time, while
the time-series panels compare the full-data fit against the dominant
68% bootstrap support.

```matlab
requestedDisplayDays = [0 2 4 6];
timeIndices = zeros(size(requestedDisplayDays));
for iDay = 1:numel(requestedDisplayDays)
    [~, timeIndices(iDay)] = min(abs(tDays - requestedDisplayDays(iDay)));
end
matchedTimeDays = tDays(timeIndices);
sigmaN = bootstrap.summary.sigma_n/f0;
sigmaS = bootstrap.summary.sigma_s/f0;
fullSigmaN = reshape(bootstrap.fullSummary.sigma_n, [], 1)/f0;
fullSigmaS = reshape(bootstrap.fullSummary.sigma_s, [], 1)/f0;
halfWidth = 1.1 * max(abs([sigmaN(:); sigmaS(:); fullSigmaN; fullSigmaS]));
thetaReference = 0.5 * unwrap(atan2(fullSigmaS, fullSigmaN));
statisticsSeries = cell(numel(tDays), 1);
strainSummarySeries = cell(numel(tDays), 1);
modeSigma = NaN(numel(tDays), 1);
sigmaBounds = NaN(numel(tDays), 2);
modeThetaDegrees = NaN(numel(tDays), 1);
thetaBoundsDegrees = NaN(numel(tDays), 2);

for iTime = 1:numel(tDays)
    samples = [reshape(sigmaN(iTime, :), [], 1), reshape(sigmaS(iTime, :), [], 1)];
    statisticsSeries{iTime} = KernelDensityEstimate.planarStatisticsFromData( ...
        samples, referencePoint=[fullSigmaN(iTime), fullSigmaS(iTime)], ...
        minimum=[-halfWidth, -halfWidth], maximum=[halfWidth, halfWidth], summaryMass=0.68);
    strainSummarySeries{iTime} = KernelDensityEstimate.strainSummaryFromPlanarStatistics( ...
        statisticsSeries{iTime}, thetaReference=thetaReference(iTime));
    modeSigma(iTime) = strainSummarySeries{iTime}.modeSigma;
    sigmaBounds(iTime, :) = strainSummarySeries{iTime}.sigmaBounds;
    modeThetaDegrees(iTime) = strainSummarySeries{iTime}.modeTheta * 180/pi;
    thetaBoundsDegrees(iTime, :) = strainSummarySeries{iTime}.thetaBounds * 180/pi;
end

fullSigma = hypot(bootstrap.fullSummary.sigma_n, bootstrap.fullSummary.sigma_s)/f0;
fullThetaDegrees = GriddedStreamfunction.visualPrincipalStrainAngle(bootstrap.fullSummary.sigma_n, bootstrap.fullSummary.sigma_s);
ringRadii = (2:2:10) * 1e-6 / f0;
contourMasses = (0.8:-0.1:0.1).';

figure(Color="w", Position=[100 100 760 860]);
tlError = tiledlayout(4, 2, TileSpacing="compact", Padding="compact");
cmap = parula(256);
cmap(1, :) = 1;
colormap(cmap);

for iPanel = 1:numel(requestedDisplayDays)
    ax = nexttile(tlError, iPanel);
    samples = [reshape(sigmaN(timeIndices(iPanel), :), [], 1), reshape(sigmaS(timeIndices(iPanel), :), [], 1)];
    KernelDensityEstimate.plotPlanarStatistics( ...
        ax, statisticsSeries{timeIndices(iPanel)}, samples=samples, contourMasses=contourMasses, ...
        showWedge=~strainSummarySeries{timeIndices(iPanel)}.containsZero, ...
        scatterColor=[0 0.4470 0.7410], wedgeColor=0.5 * [1 1 1], wedgeAlpha=0.35, ...
        wedgeResolution=100, ringRadii=ringRadii, ringColor=0.75 * [1 1 1], ...
        ringLineStyle=":", ringLineWidth=0.8);
    axis(ax, "equal")
    xlabel(ax, "\sigma_n / f_0")
    ylabel(ax, "\sigma_s / f_0")
    text(ax, 0.03, 0.97, sprintf("day %.1f", matchedTimeDays(iPanel)), Units="normalized", HorizontalAlignment="left", VerticalAlignment="top", Color=[0.2 0.2 0.2])
    box(ax, "on")
end

axSigma = nexttile(tlError, 5, [1 2]);
sigmaBand = fill( ...
    axSigma, [tDays; flipud(tDays)], [sigmaBounds(:, 2); flipud(sigmaBounds(:, 1))], ...
    [0 0.4470 0.7410], EdgeColor="none", FaceAlpha=0.22);
hold(axSigma, "on")
sigmaModeLine = plot(axSigma, tDays, modeSigma, Color=[0 0.4470 0.7410], LineWidth=1.5);
sigmaFullLine = plot(axSigma, tDays, fullSigma, "k--", LineWidth=1.3);
arrayfun(@(day) xline(axSigma, day, ":", Color=0.85 * [1 1 1], HandleVisibility="off"), matchedTimeDays);
ylabel(axSigma, "\sigma / f_0")
xlim(axSigma, [tDays(1), tDays(end)])
legend(axSigma, [sigmaModeLine, sigmaBand, sigmaFullLine], "KDE mode", "68% contour band", "full fit", Location="northwest")
box(axSigma, "on")

axTheta = nexttile(tlError, 7, [1 2]);
fill( ...
    axTheta, [tDays; flipud(tDays)], [thetaBoundsDegrees(:, 2); flipud(thetaBoundsDegrees(:, 1))], ...
    [0 0.4470 0.7410], EdgeColor="none", FaceAlpha=0.22);
hold(axTheta, "on")
plot(axTheta, tDays, modeThetaDegrees, Color=[0 0.4470 0.7410], LineWidth=1.5)
plot(axTheta, tDays, fullThetaDegrees, "k--", LineWidth=1.3)
arrayfun(@(day) xline(axTheta, day, ":", Color=0.85 * [1 1 1], HandleVisibility="off"), matchedTimeDays);
xlabel(axTheta, "time (days)")
ylabel(axTheta, "\theta (deg)")
xlim(axTheta, [tDays(1), tDays(end)])
box(axTheta, "on")

title(tlError, "Bootstrap error assessment")
```

![Whole-drifter bootstrap KDEs show how the recovered Site 1 strain state moves through `(\sigma_n,\sigma_s)` space, while the time-series panels compare the full-data fit against the dominant 68% bootstrap support.](./bootstrap-error-assessment-and-model-choice/bootstrap-error-assessment.png)

*Whole-drifter bootstrap KDEs show how the recovered Site 1 strain state moves through `(\sigma_n,\sigma_s)` space, while the time-series panels compare the full-data fit against the dominant 68% bootstrap support.*

## Compare a small set of candidate models

Model choice asks how a few nearby mesoscale model families trade
complexity against the scalar bootstrap diagnostics. Here the candidates
differ only in the mesoscale spline time dependence and in whether the
hard zero-vorticity constraint is enforced.

```matlab
candidatePsiS = [2 2 1; 2 2 3];
candidateConstraints = ["zeroVorticity", "none"];
nCandidateBootstraps = 12;
nCandidateModels = size(candidatePsiS, 1) * numel(candidateConstraints);
psiSColumn = strings(nCandidateModels, 1);
mesoscaleConstraintColumn = strings(nCandidateModels, 1);
mesoscaleDegreesOfFreedomColumn = NaN(nCandidateModels, 1);
fullFitKappaColumn = NaN(nCandidateModels, 1);
fullFitCoherenceColumn = NaN(nCandidateModels, 1);
bestFitKappaColumn = NaN(nCandidateModels, 1);
bestFitCoherenceColumn = NaN(nCandidateModels, 1);
plotLabelColumn = strings(nCandidateModels, 1);

iModel = 0;
for iConstraint = 1:numel(candidateConstraints)
    mesoscaleConstraint = candidateConstraints(iConstraint);
    for iPsiS = 1:size(candidatePsiS, 1)
        iModel = iModel + 1;
        psiS = candidatePsiS(iPsiS, :);
        candidateBootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
            trajectories, nBootstraps=nCandidateBootstraps, randomSeed=0, ...
            queryTimes=t, scoreStride=6, psiS=psiS, fastS=3, mesoscaleConstraint=mesoscaleConstraint);
        psiSColumn(iModel) = "[" + join(string(psiS), " ") + "]";
        mesoscaleConstraintColumn(iModel) = mesoscaleConstraint;
        mesoscaleDegreesOfFreedomColumn(iModel) = candidateBootstrap.mesoscaleDegreesOfFreedom;
        fullFitKappaColumn(iModel) = candidateBootstrap.fullFitKappa;
        fullFitCoherenceColumn(iModel) = candidateBootstrap.fullFitCoherence;
        bestFitKappaColumn(iModel) = candidateBootstrap.bestFitKappa;
        bestFitCoherenceColumn(iModel) = candidateBootstrap.bestFitCoherence;
        if mesoscaleConstraint == "zeroVorticity"
            plotLabelColumn(iModel) = "zeta=0 " + psiSColumn(iModel);
        else
            plotLabelColumn(iModel) = "free zeta " + psiSColumn(iModel);
        end
    end
end

plotTable = table( ...
    psiSColumn, mesoscaleConstraintColumn, mesoscaleDegreesOfFreedomColumn, ...
    fullFitKappaColumn, fullFitCoherenceColumn, bestFitKappaColumn, ...
    bestFitCoherenceColumn, plotLabelColumn, ...
    VariableNames=["psiS", "mesoscaleConstraint", "mesoscaleDegreesOfFreedom", "fullFitKappa", "fullFitCoherence", "bestFitKappa", "bestFitCoherence", "plotLabel"]);
plotTable = sortrows(plotTable, ["mesoscaleConstraint", "psiS"]);
diagnosticsTable = plotTable(:, ["psiS", "mesoscaleConstraint", "mesoscaleDegreesOfFreedom", "fullFitKappa", "fullFitCoherence", "bestFitKappa", "bestFitCoherence"]);
```

```text
      psiS       mesoscaleConstraint    mesoscaleDegreesOfFreedom    fullFitKappa    fullFitCoherence    bestFitKappa    bestFitCoherence
    _________    ___________________    _________________________    ____________    ________________    ____________    ________________

    "[2 2 1]"      "none"                          16                   1.0224           0.21936           0.65944           0.25534
    "[2 2 3]"      "none"                          32                  0.29842           0.22827           0.29154           0.22404
    "[2 2 1]"      "zeroVorticity"                  8                  0.17676           0.14883           0.13082           0.16719
    "[2 2 3]"      "zeroVorticity"                 16                  0.17433            0.1518           0.13985           0.16352
```

*The model-choice table compares the four tutorial candidates using mesoscale degrees of freedom, diffusivity, and coherence diagnostics.*

```matlab
figure(Color="w", Position=[100 100 900 360]);
tlModelChoice = tiledlayout(1, 2, TileSpacing="compact", Padding="compact");
axCoherence = nexttile;
axKappa = nexttile;
hold(axCoherence, "on")
hold(axKappa, "on")
constraintColors = struct("zeroVorticity", [0 0.4470 0.7410], "none", [0.8500 0.3250 0.0980]);

for iModel = 1:height(plotTable)
    constraint = plotTable.mesoscaleConstraint(iModel);
    if constraint == "zeroVorticity"
        marker = "o";
    else
        marker = "s";
    end
    color = constraintColors.(char(constraint));
    xValue = plotTable.mesoscaleDegreesOfFreedom(iModel);

    if isfinite(plotTable.fullFitCoherence(iModel)) && isfinite(plotTable.bestFitCoherence(iModel))
        plot(axCoherence, [xValue xValue], [plotTable.fullFitCoherence(iModel) plotTable.bestFitCoherence(iModel)], Color=0.8 * [1 1 1], LineWidth=1.2)
    end
    if isfinite(plotTable.fullFitCoherence(iModel))
        plot(axCoherence, xValue, plotTable.fullFitCoherence(iModel), marker, MarkerSize=8, MarkerFaceColor="none", MarkerEdgeColor=color, LineWidth=1.2)
    end
    if isfinite(plotTable.bestFitCoherence(iModel))
        plot(axCoherence, xValue, plotTable.bestFitCoherence(iModel), marker, MarkerSize=8, MarkerFaceColor=color, MarkerEdgeColor=color, LineWidth=1.2)
        text(axCoherence, xValue, plotTable.bestFitCoherence(iModel), "  " + plotTable.plotLabel(iModel), VerticalAlignment="middle")
    end

    if isfinite(plotTable.fullFitKappa(iModel)) && isfinite(plotTable.bestFitKappa(iModel))
        plot(axKappa, [xValue xValue], [plotTable.fullFitKappa(iModel) plotTable.bestFitKappa(iModel)], Color=0.8 * [1 1 1], LineWidth=1.2)
    end
    if isfinite(plotTable.fullFitKappa(iModel))
        plot(axKappa, xValue, plotTable.fullFitKappa(iModel), marker, MarkerSize=8, MarkerFaceColor="none", MarkerEdgeColor=color, LineWidth=1.2)
    end
    if isfinite(plotTable.bestFitKappa(iModel))
        plot(axKappa, xValue, plotTable.bestFitKappa(iModel), marker, MarkerSize=8, MarkerFaceColor=color, MarkerEdgeColor=color, LineWidth=1.2)
        text(axKappa, xValue, plotTable.bestFitKappa(iModel), "  " + plotTable.plotLabel(iModel), VerticalAlignment="middle")
    end
end

xlabel(axCoherence, "mesoscale degrees of freedom")
ylabel(axCoherence, "coherence")
box(axCoherence, "on")
title(axCoherence, "coherence")

xlabel(axKappa, "mesoscale degrees of freedom")
ylabel(axKappa, "\kappa (m^2/s)")
box(axKappa, "on")
title(axKappa, "\kappa (m^2/s)")

title(tlModelChoice, "Bootstrap model choice")
```

![A small candidate sweep compares model complexity against the scalar bootstrap diagnostics. Open markers show the full-data fit and filled markers show the top-ranked bootstrap replicate for each candidate.](./bootstrap-error-assessment-and-model-choice/bootstrap-model-choice.png)

*A small candidate sweep compares model complexity against the scalar bootstrap diagnostics. Open markers show the full-data fit and filled markers show the top-ranked bootstrap replicate for each candidate.*
