classdef NarrowEscapeProblem < KinematicModel
    % NarrowEscapeProblem Rectangular chamber with a narrow opening in the left wall.
    %
    % This model has zero deterministic velocity,
    %
    % $$ u(t,x,y) = 0, \qquad v(t,x,y) = 0, $$
    %
    % and represents the geometry through polygonal obstacles that leave a
    % gap of width `W` in the left wall. The problem is useful for testing
    % stochastic escape through a small opening.
    %
    % Topic:
    %   Models

    properties
        L = 10e3
        % Chamber length in meters.

        delta = 1e3
        % Wall thickness in meters.

        W = 3e3
        % Opening width in meters.
    end

    methods
        function self = NarrowEscapeProblem()
            % NarrowEscapeProblem Create the default narrow-escape geometry.
            self.xVisualLimits = 0.1 * (self.L + 2 * self.delta) * [-1 1] + [-self.delta self.L + self.delta];
            self.yVisualLimits = 0.1 * (self.L + 2 * self.delta) * [-1 1] + [-self.delta self.L + self.delta];
            self.visualScale = 1e3;
            self.regenerateObstacles();
            self.name = 'Narrow escape';
        end

        function set.W(self, W)
            self.W = W;
            self.regenerateObstacles();
        end

        function set.L(self, L)
            self.L = L;
            self.regenerateObstacles();
        end

        function set.delta(self, delta)
            self.delta = delta;
            self.regenerateObstacles();
        end
    end
end
