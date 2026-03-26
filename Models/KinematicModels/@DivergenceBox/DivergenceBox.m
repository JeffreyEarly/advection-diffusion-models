classdef DivergenceBox < KinematicModel
    % DivergenceBox Box with alternating Gaussian convergence and divergence cells.
    %
    % The model defines a steady velocity field by summing Gaussian cells
    % centered at `(x_c, y_c)`:
    %
    % $$ u(x,y) = \sum_{i,j} s_i e^{1/2} U \frac{x - x_c}{L_g} \exp\left(-\frac{r^2}{2L_g^2}\right), $$
    %
    % $$ v(x,y) = \sum_{i,j} s_i e^{1/2} U \frac{y - y_c}{L_g} \exp\left(-\frac{r^2}{2L_g^2}\right), $$
    %
    % where `s_i = (-1)^i` alternates the sign of neighboring cells.
    %
    % Topic:
    %   Models

    properties
        Lx = 10e3
        % Domain width in meters.

        Ly = 5e3
        % Domain height in meters.

        Lg = 0.5e3
        % Gaussian length scale in meters.

        n = 2
        % Number of cells along x.

        m = 1
        % Number of cells along y.

        U = 1.0
        % Velocity scale in m s^-1.
    end

    methods
        function self = DivergenceBox()
            % DivergenceBox Create the default divergence-box model.
            self.xlim = [0 self.Lx];
            self.ylim = [0 self.Ly];

            epsilonX = 0.05 * self.Lx;
            epsilonY = 0.05 * self.Ly;
            self.xVisualLimits = [0 - epsilonX self.Lx + epsilonX];
            self.yVisualLimits = [0 - epsilonY self.Ly + epsilonY];
            self.visualScale = 1e3;
            self.name = 'Divergence box';
        end

        function uValue = u(self, t, x, y)
            % u Evaluate the x-velocity component.
            uValue = zeros(size(x));
            for i = 1:self.n
                for j = 1:self.m
                    xc = (2 * i - 1) * self.Lx / (self.n + 2);
                    yc = j * self.Ly / (self.m + 1);
                    r2 = (x - xc).^2 + (y - yc).^2;
                    gaussian = exp(-r2 / (2 * self.Lg * self.Lg));
                    uValue = uValue + ((-1).^i) * exp(0.5) * self.U * ((x - xc) / self.Lg) .* gaussian;
                end
            end
        end

        function vValue = v(self, t, x, y)
            % v Evaluate the y-velocity component.
            vValue = zeros(size(x));
            for i = 1:self.n
                for j = 1:self.m
                    xc = (2 * i - 1) * self.Lx / (self.n + 2);
                    yc = j * self.Ly / (self.m + 1);
                    r2 = (x - xc).^2 + (y - yc).^2;
                    gaussian = exp(-r2 / (2 * self.Lg * self.Lg));
                    vValue = vValue + ((-1).^i) * exp(0.5) * self.U * ((y - yc) / self.Lg) .* gaussian;
                end
            end
        end
    end
end
