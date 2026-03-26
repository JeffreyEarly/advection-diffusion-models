function varargout = plotStreamfunction(self, t)
% plotStreamfunction Plot streamfunction contours on the current axes.
%
% Declaration:
%   `plotStreamfunction(self)`
%   `plotStreamfunction(self, t)`
%
% Parameters:
%   `t` - Scalar plotting time in seconds. The default is `0`.
%
% Returns:
%   `X`, `Y` - Optional plotting grids in meters.
arguments
    self (1,1) StreamfunctionModel
    t (1,1) double = 0
end

numPoints = 150;
xg = linspace(min(self.xVisualLimits), max(self.xVisualLimits), 2 * numPoints).';
yg = linspace(min(self.yVisualLimits), max(self.yVisualLimits), numPoints).';
[Xf, Yf] = meshgrid(xg, yg);

if ~isempty(self.obstacles)
    x = reshape(Xf, [], 1);
    y = reshape(Yf, [], 1);
    mask = ~isinterior(self.obstacles, x, y);
    mask = reshape(mask, size(Xf));
    psiValues = self.psi(t, Xf(mask), Yf(mask));
    levels = linspace(min(psiValues), max(psiValues), 20);

    contour(Xf / self.visualScale, Yf / self.visualScale, mask .* self.psi(t, Xf, Yf), levels)
    hold on
    plot(scale(self.obstacles, 1 / self.visualScale))
else
    contour(Xf / self.visualScale, Yf / self.visualScale, self.psi(t, Xf, Yf))
end

axis equal
xlim([min(self.xVisualLimits) max(self.xVisualLimits)] / self.visualScale)
ylim([min(self.yVisualLimits) max(self.yVisualLimits)] / self.visualScale)
xlabel('km')
ylabel('km')
title(sprintf('%s', self.name))

if nargout == 2
    varargout{1} = Xf;
    varargout{2} = Yf;
end
