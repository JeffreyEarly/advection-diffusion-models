scriptDir = fileparts(mfilename('fullpath'));
siteFit = loadGriddedStreamfunctionSiteFit(1);

scaleFactor = 1;
run(fullfile(scriptDir, 'LoadFigureDefaults.m'));

tDays = siteFit.t/86400;

figure('Units', 'points', 'Position', [50 50 figureWidthMedium 420*scaleFactor])
set(gcf, 'PaperPositionMode', 'auto')
set(gcf, 'Color', 'w');

tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
plot(tDays, siteFit.sigma/siteFit.f0, 'LineWidth', 2)
ylabel('\sigma / f_0', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
xlim([tDays(1) tDays(end)])
set(gca, 'FontSize', figureAxisTickSize, 'FontName', figureFont)
title('Site 1 Gridded Streamfunction Fit', 'FontSize', figureTitleSize, 'FontName', figureFont)

nexttile
plot(tDays, siteFit.thetaDegrees, 'LineWidth', 2)
ylabel('\theta (°)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
xlim([tDays(1) tDays(end)])
set(gca, 'FontSize', figureAxisTickSize, 'FontName', figureFont)

nexttile
plot(tDays, siteFit.uBackground, 'LineWidth', 2)
hold on
plot(tDays, siteFit.vBackground, 'LineWidth', 2)
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
plot(siteFit.qObserved/1000, siteFit.rObserved/1000, 'Color', 0.8*[1 1 1])
hold on
plot(siteFit.qMesoscale/1000, siteFit.rMesoscale/1000, 'LineWidth', 1.5)
axis equal
xlabel('q (km)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
ylabel('r (km)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
title('Centered trajectories', 'FontSize', figureTitleSize, 'FontName', figureFont)
legend({'Observed', 'Mesoscale fit'}, 'Location', 'best', 'FontSize', figureLegendSize)
set(gca, 'FontSize', figureAxisTickSize, 'FontName', figureFont)

nexttile
plot(tDays, siteFit.observedSpeedRms, 'Color', 0.2*[1 1 1], 'LineWidth', 1.5)
hold on
plot(tDays, siteFit.mesoscaleSpeedRms, 'LineWidth', 2)
plot(tDays, siteFit.backgroundSpeedRms, 'LineWidth', 2)
plot(tDays, siteFit.submesoscaleSpeedRms, 'LineWidth', 2)
xlabel('time (days)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
ylabel('RMS speed (m/s)', 'FontSize', figureAxisLabelSize, 'FontName', figureFont)
xlim([tDays(1) tDays(end)])
legend({'Observed', 'Mesoscale', 'Background', 'Submesoscale'}, 'Location', 'best', 'FontSize', figureLegendSize)
set(gca, 'FontSize', figureAxisTickSize, 'FontName', figureFont)
