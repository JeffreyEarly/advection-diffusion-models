scriptDir = locateExampleScriptDirectory();
timeDays = [0.5 2.5 5.0];
timeDays = [3 4 5];
nBootstraps = 1000;
randomSeed = 0;
shouldSaveFigures = 0;
scaleFactor = 1;

requireMethodComparisonFunction("buildSite1OscroftFairComparison")
requireMethodComparisonFunction("buildSite1BootstrapStrainKdeComparison")

comparison = buildSite1OscroftFairComparison( ...
    nBootstraps=nBootstraps, ...
    randomSeed=randomSeed, ...
    shouldLoadPublishedBaseline=false, ...
    displaySummary=false);
kdeComparison = buildSite1BootstrapStrainKdeComparison(comparison, timeDays);

run(fullfile(scriptDir, "LoadFigureDefaults.m"));

disp(table( ...
    kdeComparison.requestedTimeDays, ...
    kdeComparison.matchedTimeDays, ...
    kdeComparison.timeIndices, ...
    VariableNames=["RequestedDay", "MatchedDay", "QueryIndex"]))

methodColors = [0.8500 0.3250 0.0980; 0 0.4470 0.7410];
figureHeight = max(360, 170 * numel(kdeComparison.requestedTimeDays)) * scaleFactor;

figure(Color="w", Units="points", Position=[50 50 figureWidthMedium figureHeight]);
tl = tiledlayout(numel(kdeComparison.requestedTimeDays), 2, TileSpacing="compact", Padding="compact");
cmap = parula(256);
cmap(1, :) = 1;
colormap(cmap);

for iTime = 1:numel(kdeComparison.requestedTimeDays)
    for iMethod = 1:numel(kdeComparison.methodLabels)
        ax = nexttile(tl);
        panel = kdeComparison.panels(iTime, iMethod);
        plotStrainKdePanel(ax, panel, kdeComparison.bounds, methodColors(iMethod, :));

        xlabel(ax, "\sigma_n (s^{-1})", FontSize=figureAxisLabelSize, FontName=figureFont);
        ylabel(ax, "\sigma_s (s^{-1})", FontSize=figureAxisLabelSize, FontName=figureFont);
        set(ax, FontSize=figureAxisTickSize, FontName=figureFont);

        if iTime == 1
            title(ax, kdeComparison.methodLabels(iMethod), FontSize=figureTitleSize, FontName=figureFont);
        end

        text(ax, 0.03, 0.97, sprintf("t = %.3f d\nindex %d", panel.matchedTimeDays, panel.timeIndex), ...
            Units="normalized", HorizontalAlignment="left", VerticalAlignment="top", ...
            FontSize=figureLegendSize, FontName=figureFont, Color=[0.2 0.2 0.2]);
    end
end

title(tl, "Site 1 fair-comparison bootstrap strain KDE", FontSize=figureTitleSize, FontName=figureFont);

if shouldSaveFigures == 1
    print(fullfile(scriptDir, "Site1OscroftFairComparisonStrainKDE.eps"), "-depsc2");
end

function plotStrainKdePanel(ax, panel, bounds, panelColor)
hold(ax, "on");
drawReferenceRings(ax);
contourf(ax, panel.gridVectors{1}, panel.gridVectors{2}, panel.density.', panel.contourLevels, LineStyle="none");
scatter(ax, panel.samples(:, 1), panel.samples(:, 2), 12, panelColor, "filled", ...
    MarkerEdgeColor="none", MarkerFaceAlpha=0.18);
scatter(ax, panel.fullPoint(1), panel.fullPoint(2), 42, "w", "filled", MarkerEdgeColor="k", LineWidth=1.1);
axis(ax, "equal");
xlim(ax, bounds.minimum(1) + [0, 1] * (bounds.maximum(1) - bounds.minimum(1)));
ylim(ax, bounds.minimum(2) + [0, 1] * (bounds.maximum(2) - bounds.minimum(2)));
box(ax, "on");
end

function drawReferenceRings(ax)
for radius = (2:2:10) * 1e-6
    rectangle(ax, Position=[-radius -radius 2 * radius 2 * radius], Curvature=[1 1], ...
        EdgeColor=0.75 * [1 1 1], LineStyle=":", LineWidth=0.8);
end
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
