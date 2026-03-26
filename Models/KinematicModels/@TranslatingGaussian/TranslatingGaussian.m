classdef TranslatingGaussian < StreamfunctionModel
    % TranslatingGaussian Translating Gaussian eddy streamfunction.
    %
    % The streamfunction is
    %
    % $$ \psi(t,x,y) = e^{1/2} U L \exp\left(-\frac{(x - c_x t)^2 + (y - c_y t)^2}{2L^2}\right), $$
    %
    % which produces a translating Gaussian vortex with peak velocity scale
    % `U` and translation velocity `(cx, cy)`.
    %
    % Topic:
    %   Models

    properties
        L = 60e3
        % Eddy length scale in meters.

        U = 0.12
        % Peak velocity scale in m s^-1.

        cx = -0.0267
        % Translation speed in x in m s^-1.

        cy = 0
        % Translation speed in y in m s^-1.
    end

    methods
        function self = TranslatingGaussian()
            % TranslatingGaussian Create the default translating Gaussian model.
            self.xVisualLimits = 3 * self.L * [-1 1];
            self.yVisualLimits = 3 * self.L * [-1 1];
            self.visualScale = 1e3;
            self.name = 'Translating Gaussian';
        end

        function psiValue = psi(self, t, x, y)
            % psi Evaluate the Gaussian streamfunction.
            r2 = (x - self.cx * t).^2 + (y - self.cy * t).^2;
            gaussian = exp(-r2 / (2 * self.L * self.L));
            psiValue = exp(0.5) * self.U * self.L * gaussian;
        end

        function uValue = u(self, t, x, y)
            % u Evaluate the x-velocity component.
            r2 = (x - self.cx * t).^2 + (y - self.cy * t).^2;
            gaussian = exp(-r2 / (2 * self.L * self.L));
            uValue = exp(0.5) * self.U * ((y - self.cy * t) / self.L) .* gaussian;
        end

        function vValue = v(self, t, x, y)
            % v Evaluate the y-velocity component.
            r2 = (x - self.cx * t).^2 + (y - self.cy * t).^2;
            gaussian = exp(-r2 / (2 * self.L * self.L));
            vValue = -exp(0.5) * self.U * ((x - self.cx * t) / self.L) .* gaussian;
        end
    end
end
