scriptDir = fileparts(mfilename("fullpath"));
siteFit = loadGriddedStreamfunctionSiteFit(1);

scaleFactor = 1;
run(fullfile(scriptDir, "LoadFigureDefaults.m"));

tDays = siteFit.t/86400;
xLimits = [tDays(1), tDays(end)];
sigmaFixed = 0.0591 * siteFit.f0;
thetaFixed = -27.8;
thetaPlotDegrees = principalThetaDegrees(siteFit.thetaDegrees);
figureHeight = 367;

figure("Units", "points", "Position", [50 50 figureWidthMedium 450*scaleFactor])
set(gcf, "PaperPositionMode", "auto")
set(gcf, "Color", "w");

tl = tiledlayout(4, 1, "TileSpacing", "compact", "Padding", "compact");

ax1 = nexttile;
hold(ax1, "on")
plot(ax1, xLimits, sigmaFixed*[1 1]/siteFit.f0, "k--", "LineWidth", 1.0*scaleFactor)
plot(ax1, tDays, siteFit.sigma/siteFit.f0, "LineWidth", 2)
ylabel(ax1, "\sigma (f_0)", "FontSize", figureAxisLabelSize, "FontName", figureFont)
xlim(ax1, xLimits)
styleAxis(ax1, figureFont, figureAxisTickSize)

ax2 = nexttile;
hold(ax2, "on")
plot(ax2, xLimits, thetaFixed*[1 1], "k--", "LineWidth", 1.0*scaleFactor)
plot(ax2, tDays, thetaPlotDegrees, "LineWidth", 2)
ylabel(ax2, "\theta (°)", "FontSize", figureAxisLabelSize, "FontName", figureFont)
xlim(ax2, xLimits)
ylim(ax2, [-100 100])
yticks(ax2, [-90 -45 0 45 90])
styleAxis(ax2, figureFont, figureAxisTickSize)

ax3 = nexttile;
hold(ax3, "on")
u0Plot = plot(ax3, tDays, siteFit.effective.uBackground, "LineWidth", 2);
v0Plot = plot(ax3, tDays, siteFit.effective.vBackground, "LineWidth", 2);
ylabel(ax3, "u_0/v_0 (m/s)", "FontSize", figureAxisLabelSize, "FontName", figureFont)
xlim(ax3, xLimits)
ylim(ax3, [-0.25 0.25])
legend(ax3, [u0Plot, v0Plot], "u_0", "v_0", "Location", "southwest", ...
    "NumColumns", 2, "FontSize", figureLegendSize, "Box", "off")
styleAxis(ax3, figureFont, figureAxisTickSize)

ax4 = nexttile;
basisColors = gray(size(siteFit.timeBasisMatrix, 2) + 4);
ax4.ColorOrder = basisColors(3:(size(siteFit.timeBasisMatrix, 2) + 2), :);
plot(ax4, tDays, siteFit.timeBasisMatrix, "LineWidth", 2)
xlim(ax4, xLimits)
ax4.YTick = [];
xlabel(ax4, "time (days)", "FontSize", figureAxisLabelSize, "FontName", figureFont)
text(ax4, 0.02, 0.82, "Matched cubic dof = 4 basis", ...
    "Units", "normalized", "FontSize", figureAxisLabelSize, "FontName", figureFont)
text(ax4, 0.02, 0.62, "Bootstrap error bars omitted in gridded example", ...
    "Units", "normalized", "FontSize", figureAxisTickSize, "FontName", figureFont)
styleAxis(ax4, figureFont, figureAxisTickSize)

title(tl, "Site 1 gridded streamfunction fit")

figure("Units", "points", "Position", [50 50 figureWidthOneColumn figureHeight*scaleFactor])
set(gcf, "PaperPositionMode", "auto")
set(gcf, "Color", "w");

axMeso = axes;
plot(axMeso, siteFit.effective.xMesoscale/1000, siteFit.effective.yMesoscale/1000, "LineWidth", 1.5)
axis(axMeso, "equal")
xlabel(axMeso, "km", "FontSize", figureAxisLabelSize, "FontName", figureFont)
ylabel(axMeso, "km", "FontSize", figureAxisLabelSize, "FontName", figureFont)
text(axMeso, 0.05, 0.92, "meso (effective)", ...
    "Units", "normalized", "FontSize", figureAxisLabelSize, "FontName", figureFont)
styleAxis(axMeso, figureFont, figureAxisTickSize)

figure("Units", "points", "Position", [50 50 figureWidthOneColumn figureHeight*scaleFactor])
set(gcf, "PaperPositionMode", "auto")
set(gcf, "Color", "w");

tlDecomp = tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact");

axBackground = nexttile;
plot(axBackground, siteFit.effective.xBackground/1000, siteFit.effective.yBackground/1000, "LineWidth", 1.5)
axis(axBackground, "equal")
axBackground.YAxisLocation = "right";
axBackground.XTickLabel = [];
ylabel(axBackground, "km", "FontSize", figureAxisLabelSize, "FontName", figureFont)
text(axBackground, 0.05, 0.88, "bg", ...
    "Units", "normalized", "FontSize", figureAxisLabelSize, "FontName", figureFont)
styleAxis(axBackground, figureFont, figureAxisTickSize)

axSubmeso = nexttile;
plot(axSubmeso, siteFit.effective.xSubmesoscale/1000, siteFit.effective.ySubmesoscale/1000, "LineWidth", 1.5)
axis(axSubmeso, "equal")
axSubmeso.YAxisLocation = "right";
xlabel(axSubmeso, "km", "FontSize", figureAxisLabelSize, "FontName", figureFont)
ylabel(axSubmeso, "km", "FontSize", figureAxisLabelSize, "FontName", figureFont)
text(axSubmeso, 0.05, 0.88, "sm", ...
    "Units", "normalized", "FontSize", figureAxisLabelSize, "FontName", figureFont)
styleAxis(axSubmeso, figureFont, figureAxisTickSize)

[decompositionXLimits, decompositionYLimits] = sharedTrajectoryLimits( ...
    {siteFit.effective.xBackground, siteFit.effective.xSubmesoscale}, ...
    {siteFit.effective.yBackground, siteFit.effective.ySubmesoscale}, 1000);
xlim(axBackground, decompositionXLimits)
xlim(axSubmeso, decompositionXLimits)
ylim(axBackground, decompositionYLimits)
ylim(axSubmeso, decompositionYLimits)

title(tlDecomp, "Site 1 background and submesoscale decomposition")

figure("Units", "points", "Position", [50 50 figureWidthTwoColumn figureHeight*scaleFactor])
set(gcf, "PaperPositionMode", "auto")
set(gcf, "Color", "w");

axCom = axes;
plot(axCom, siteFit.effective.qMesoscale/1000, siteFit.effective.rMesoscale/1000, "LineWidth", 1.5)
axis(axCom, "equal")
xlabel(axCom, "km", "FontSize", figureAxisLabelSize, "FontName", figureFont)
ylabel(axCom, "km", "FontSize", figureAxisLabelSize, "FontName", figureFont)
text(axCom, 0.05, 0.92, "meso (centre-of-mass, effective)", ...
    "Units", "normalized", "FontSize", figureAxisLabelSize, "FontName", figureFont)
styleAxis(axCom, figureFont, figureAxisTickSize)

function styleAxis(ax, figureFont, figureAxisTickSize)
ax.FontSize = figureAxisTickSize;
ax.FontName = figureFont;
box(ax, "on")
end

function [xLimits, yLimits] = sharedTrajectoryLimits(xCell, yCell, scale)
xValues = cellfun(@(values) values(:), xCell, 'UniformOutput', false);
yValues = cellfun(@(values) values(:), yCell, 'UniformOutput', false);
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
limits = limits + 0.05*range*[-1 1];
end

function thetaDegrees = principalThetaDegrees(thetaDegrees)
thetaDegrees = thetaDegrees(:);
thetaDegrees = mod(thetaDegrees + 45, 90) - 45;
end
