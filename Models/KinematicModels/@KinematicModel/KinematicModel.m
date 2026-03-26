classdef KinematicModel < handle
    % KinematicModel Base class for two-dimensional deterministic velocity fields.
    %
    % A kinematic model defines particle motion through
    %
    % $$ \frac{dx}{dt} = u(t,x,y), \qquad \frac{dy}{dt} = v(t,x,y). $$
    %
    % Subclasses specify the velocity field itself and, optionally, finite
    % domain limits, periodic directions, polygonal obstacles, and plotting
    % limits. The base class provides common geometry filtering and
    % visualization helpers.
    %
    % Topic:
    %   Models

    properties
        xlim (1,2) double = [-Inf Inf]
        % Finite or infinite x-domain limits in meters.

        ylim (1,2) double = [-Inf Inf]
        % Finite or infinite y-domain limits in meters.

        isXPeriodic (1,1) logical = false
        % True when the x-direction wraps periodically across `xlim`.

        isYPeriodic (1,1) logical = false
        % True when the y-direction wraps periodically across `ylim`.

        obstacles = []
        % Polygonal obstacle geometry as a `polyshape` array.

        isSpherical (1,1) logical = false
        % True when the x-coordinate represents longitude on a sphere.

        xVisualLimits (1,2) double = [-1 1]
        % Finite x-limits in meters used for plotting.

        yVisualLimits (1,2) double = [-1 1]
        % Finite y-limits in meters used for plotting.

        visualScale (1,1) double {mustBePositive} = 1e3
        % Plotting scale factor. A value of `1e3` renders meters as km.

        name = ''
        % Descriptive model name used in figures and docs.
    end

    properties (Access = private)
        uFunction = []
        vFunction = []
    end

    methods
        function uValue = u(self, t, x, y)
            % u Evaluate the x-velocity component.
            %
            % Declaration:
            %   `uValue = u(self, t, x, y)`
            %
            % Parameters:
            %   `t` - Scalar evaluation time in seconds.
            %
            %   `x`, `y` - Position arrays in meters. Output shape matches
            %   the input shape.
            %
            % Returns:
            %   `uValue` - x-velocity in m s^-1.
            if isempty(self.uFunction)
                uValue = zeros(size(x));
            else
                uValue = self.uFunction(t, x, y);
            end
        end

        function vValue = v(self, t, x, y)
            % v Evaluate the y-velocity component.
            %
            % Declaration:
            %   `vValue = v(self, t, x, y)`
            %
            % Parameters and returns follow `u(self, t, x, y)`.
            if isempty(self.vFunction)
                vValue = zeros(size(x));
            else
                vValue = self.vFunction(t, x, y);
            end
        end
    end

    methods (Static)
        function bounds = boundsFromLimits(xlim, ylim)
            % boundsFromLimits Convert rectangular limits to polygon vertices.
            %
            % Declaration:
            %   `bounds = boundsFromLimits(xlim, ylim)`
            %
            % Parameters:
            %   `xlim`, `ylim` - Two-element vectors defining rectangle
            %   limits in meters.
            %
            % Returns:
            %   `bounds` - Structure with fields `xv` and `yv`.
            arguments
                xlim (1,2) double
                ylim (1,2) double
            end

            xv = [min(xlim) min(xlim) max(xlim) max(xlim)];
            yv = [max(ylim) min(ylim) min(ylim) max(ylim)];
            bounds = struct('xv', xv, 'yv', yv);
        end
    end
end
