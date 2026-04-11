scriptDir = locateExampleScriptDirectory();
shouldSaveFigures = 0;
scaleFactor = 1;

requireMethodComparisonFunction("buildSite1OscroftFairComparison")

comparison = buildSite1OscroftFairComparison();
run(fullfile(scriptDir, "LoadFigureDefaults.m"));

tDays = comparison.site.queryTimes / 86400;
f0 = comparison.site.f0;
oldFair = comparison.old.shared;
newFair = comparison.new.shared;

oldColor = [0.8500 0.3250 0.0980];
newColor = [0 0.4470 0.7410];
lineWidth = 1.5;

figure(Color="w", Units="points", Position=[50 50 figureWidthMedium 420]);
tl = tiledlayout(3, 1, TileSpacing="compact", Padding="compact");

axSigma = nexttile;
plotComparisonSeries(axSigma, tDays, hypot(oldFair.commonSummary.sigma_n, oldFair.commonSummary.sigma_s)/f0, oldFair.scores.common, oldColor, "-");
plotComparisonSeries(axSigma, tDays, hypot(newFair.commonSummary.sigma_n, newFair.commonSummary.sigma_s)/f0, newFair.scores.common, newColor, "-");
ylabel(axSigma, "\sigma (f_0)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(axSigma, [tDays(1), tDays(end)]);
set(axSigma, FontSize=figureAxisTickSize, FontName=figureFont);
box(axSigma, "on");
title(axSigma, "Site 1 fair bootstrap comparison", FontSize=figureTitleSize, FontName=figureFont);

axTheta = nexttile;
plotComparisonSeries(axTheta, tDays, oldFair.commonSeries.theta, oldFair.scores.common, oldColor, "-");
plotComparisonSeries(axTheta, tDays, newFair.commonSeries.theta, newFair.scores.common, newColor, "-");
ylabel(axTheta, "\theta (deg)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(axTheta, [tDays(1), tDays(end)]);
ylim(axTheta, [-100 100]);
yticks(axTheta, [-90 -45 0 45 90]);
set(axTheta, FontSize=figureAxisTickSize, FontName=figureFont);
box(axTheta, "on");

axVelocity = nexttile;
hold(axVelocity, "on");
oldU = plotComparisonSeries(axVelocity, tDays, oldFair.commonSummary.uCenter, oldFair.scores.common, oldColor, "-");
oldV = plotComparisonSeries(axVelocity, tDays, oldFair.commonSummary.vCenter, oldFair.scores.common, oldColor, "--");
newU = plotComparisonSeries(axVelocity, tDays, newFair.summary.uCenter, newFair.scores.common, newColor, "-");
newV = plotComparisonSeries(axVelocity, tDays, newFair.summary.vCenter, newFair.scores.common, newColor, "--");
xlabel(axVelocity, "time (days)", FontSize=figureAxisLabelSize, FontName=figureFont);
ylabel(axVelocity, "u_c/v_c (m/s)", FontSize=figureAxisLabelSize, FontName=figureFont);
xlim(axVelocity, [tDays(1), tDays(end)]);
ylim(axVelocity, [-0.25 0.25]);
legend(axVelocity, [oldU, oldV, newU, newV], ...
    "Oscroft u_c", "Oscroft v_c", "Gridded u_c", "Gridded v_c", ...
    Location="southwest", FontSize=figureLegendSize, NumColumns=2);
set(axVelocity, FontSize=figureAxisTickSize, FontName=figureFont);
box(axVelocity, "on");

title(tl, "Site 1 shared-bootstrap sanity check", FontSize=figureTitleSize, FontName=figureFont);

if shouldSaveFigures == 1
    print(fullfile(scriptDir, "Site1OscroftFairComparison.eps"), "-depsc2");
end

function medianLine = plotComparisonSeries(ax, tDays, values, scores, baseColor, lineStyle)
rankFractions = [0.90 0.68];
rankAlphas = [0.08 0.14];
quantileAlpha = 0.22;

hold(ax, "on");
for iRank = 1:numel(rankFractions)
    indices = rankedIndices(scores, rankFractions(iRank));
    envelope = rankedEnvelope(values(:, indices));
    fill(ax, [tDays; flipud(tDays)], envelope, baseColor, ...
        EdgeColor="none", FaceAlpha=rankAlphas(iRank), HandleVisibility="off");
end

quantiles = quantile(values, [0.16 0.5 0.84], 2);
fill(ax, [tDays; flipud(tDays)], [quantiles(:, 3); flipud(quantiles(:, 1))], baseColor, ...
    EdgeColor="none", FaceAlpha=quantileAlpha, HandleVisibility="off");
medianLine = plot(ax, tDays, quantiles(:, 2), Color=baseColor, LineStyle=lineStyle, LineWidth=1.5);
end

function indices = rankedIndices(scores, fraction)
[~, rankedIndices] = sort(scores, "descend");
indices = rankedIndices(1:max(1, round(fraction * numel(scores))));
end

function envelope = rankedEnvelope(values)
envelope = [max(values, [], 2); flipud(min(values, [], 2))];
end

function scriptDir = locateExampleScriptDirectory()
scriptPath = string(mfilename("fullpath"));
if strlength(scriptPath) == 0
    scriptName = string(mfilename);
    if strlength(scriptName) > 0
        scriptPath = string(which(scriptName));
    end
end

if strlength(scriptPath) == 0 && usejava("desktop")
    scriptPath = string(matlab.desktop.editor.getActiveFilename);
end

if strlength(scriptPath) == 0
    error("GriddedStreamfunction:UnableToLocateExampleScript", ...
        "Could not determine the example script location. Run the saved script file or add Estimation/MethodComparison to the MATLAB path first.");
end

scriptDir = char(fileparts(scriptPath));
end

function requireMethodComparisonFunction(functionName)
if exist(functionName, "file") == 2
    return
end

error("GriddedStreamfunction:MissingMethodComparisonDependency", ...
    "%s is not on the MATLAB path. Reload the AdvectionDiffusionModels package path so Estimation/MethodComparison is available.", functionName);
end
