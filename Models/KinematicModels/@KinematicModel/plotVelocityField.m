function varargout = plotVelocityField(self, options)
% plotVelocityField Plot the model velocity field with quivers.
%
% Declaration:
%   `plotVelocityField(self)`
%   `plotVelocityField(self, t=..., quiverScale=..., numPoints=...)`
%
% Parameters:
%   `t` - Scalar time in seconds.
%
%   `quiverScale` - Scale factor forwarded to `quiver`.
%
%   `numPoints` - Number of points in the y-direction. The x-direction
%   uses `2*numPoints` to preserve the typical wide aspect ratio.
%
% Returns:
%   `X`, `Y` - Optional plotting grids in meters.
arguments
    self (1,1) KinematicModel
    options.t (1,1) double = 0
    options.quiverScale (1,1) double = 1
    options.numPoints (1,1) double {mustBePositive, mustBeInteger} = 150
end

xg = linspace(min(self.xVisualLimits), max(self.xVisualLimits), 2 * options.numPoints).';
yg = linspace(min(self.yVisualLimits), max(self.yVisualLimits), options.numPoints).';
[X, Y] = meshgrid(xg, yg);

if ~isempty(self.obstacles)
    x = reshape(X, [], 1);
    y = reshape(Y, [], 1);
    mask = true(size(x));
    for iObstacle = 1:length(self.obstacles)
        mask = mask & ~isinterior(self.obstacles(iObstacle), x, y);
    end
    mask = reshape(mask, size(X));

    quiver(X / self.visualScale, Y / self.visualScale, mask .* self.u(options.t, X, Y), mask .* self.v(options.t, X, Y), options.quiverScale);
    hold on
    plot(scale(self.obstacles, 1 / self.visualScale))
else
    quiver(X / self.visualScale, Y / self.visualScale, self.u(options.t, X, Y), self.v(options.t, X, Y), options.quiverScale)
end

axis equal
xlim([min(self.xVisualLimits) max(self.xVisualLimits)] / self.visualScale)
ylim([min(self.yVisualLimits) max(self.yVisualLimits)] / self.visualScale)
xlabel('km')
ylabel('km')

if nargout == 2
    varargout{1} = X;
    varargout{2} = Y;
end
