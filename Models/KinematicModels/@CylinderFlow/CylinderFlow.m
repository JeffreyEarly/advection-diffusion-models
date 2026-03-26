classdef CylinderFlow < StreamfunctionModel
    % CylinderFlow Potential flow around a circular cylinder.
    %
    % The model streamfunction is
    %
    % $$ \psi(t,x,y) = U y \left(1 - \frac{R^2}{x^2 + y^2}\right), $$
    %
    % which gives the classical incompressible potential flow past a
    % cylinder of radius `R`. The cylinder interior is represented as a
    % polygonal obstacle.
    %
    % Topic:
    %   Models

    properties
        R = 60e3
        % Cylinder radius in meters.

        U = 100 * (1e3 / 86400)
        % Far-field speed in m s^-1.
    end

    methods
        function self = CylinderFlow()
            % CylinderFlow Create the default cylinder-flow model.
            self.xlim = 3 * self.R * [-1 1];
            self.ylim = 3 * self.R * [-1 1];
            self.isXPeriodic = true;

            self.xVisualLimits = 3 * self.R * [-1 1];
            self.yVisualLimits = 3 * self.R * [-1 1];
            self.visualScale = 1e3;

            theta = linspace(0, 2 * pi - 2 * pi / 30, 30);
            self.obstacles = polyshape(self.R * cos(theta), self.R * sin(theta));
            self.name = 'Cylinder flow';
        end

        function psiValue = psi(self, t, x, y)
            % psi Evaluate the potential-flow streamfunction.
            r2 = x.^2 + y.^2;
            psiValue = self.U * y .* (1 - self.R * self.R ./ r2);
        end

        function uValue = u(self, t, x, y)
            % u Evaluate the x-velocity component.
            r2 = x.^2 + y.^2;
            uValue = self.U - self.U * self.R * self.R * (x.^2 - y.^2) ./ (r2 .* r2);
        end

        function vValue = v(self, t, x, y)
            % v Evaluate the y-velocity component.
            r2 = x.^2 + y.^2;
            vValue = -self.U * self.R * self.R * (2 * x .* y) ./ (r2 .* r2);
        end
    end
end
