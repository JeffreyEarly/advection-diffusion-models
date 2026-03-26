classdef SimpleBox < KinematicModel
    % SimpleBox Rectangular box with zero deterministic velocity.
    %
    % This model defines
    %
    % $$ u(t,x,y) = 0, \qquad v(t,x,y) = 0, $$
    %
    % on a finite rectangular domain. Any motion comes from the external
    % integrator, typically through the diffusivity term.
    %
    % Topic:
    %   Models

    properties
        Lx = 10e3
        % Domain width in meters.

        Ly = 5e3
        % Domain height in meters.
    end

    methods
        function self = SimpleBox()
            % SimpleBox Create the default box model.
            self.xlim = [0 self.Lx];
            self.ylim = [0 self.Ly];

            epsilonX = 0.05 * self.Lx;
            epsilonY = 0.05 * self.Ly;
            self.xVisualLimits = [0 - epsilonX self.Lx + epsilonX];
            self.yVisualLimits = [0 - epsilonY self.Ly + epsilonY];
            self.visualScale = 1e3;
            self.name = 'No-flow box';
        end

        function uValue = u(self, t, x, y)
            % u Evaluate the zero x-velocity field.
            uValue = zeros(size(x));
        end

        function vValue = v(self, t, x, y)
            % v Evaluate the zero y-velocity field.
            vValue = zeros(size(y));
        end
    end
end
