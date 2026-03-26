classdef MeanderingJet < StreamfunctionModel
    % MeanderingJet Kinematic streamfunction model for a meandering jet.
    %
    % This model uses the Bower meandering-jet streamfunction
    %
    % $$ \psi(t,x,y) = UL\left[1 - \tanh(\gamma)\right], $$
    %
    % with
    %
    % $$ \gamma = \frac{y - A\cos(\theta)}{L\sqrt{1 + (kA\sin(\theta))^2}}, \qquad \theta = k(x - c_x t), \qquad k = \frac{2\pi}{L_x}. $$
    %
    % The resulting velocity field is incompressible and periodic in the
    % x-direction over the interval `[0 2*Lx]`.
    %
    % References:
    %   Amy Bower, "A simple kinematic mechanism for mixing fluid parcels
    %   across a meandering jet."
    %
    % Topic:
    %   Models

    properties
        L = 40e3
        % Jet half-width in meters.

        U = 100 * (1e3 / 86400)
        % Velocity scale in m s^-1.

        A = 50e3
        % Meander amplitude in meters.

        Lx = 350e3
        % Meander wavelength in meters.

        cx = 0 * 10 * (1e3 / 86400)
        % Meander phase speed in m s^-1.
    end

    properties (Dependent)
        k
        % Meander wavenumber in m^-1.
    end

    methods
        function self = MeanderingJet()
            % MeanderingJet Create the default meandering jet model.
            self.xlim = [0 2 * self.Lx];
            self.ylim = [-5 * self.L 5 * self.L];
            self.isXPeriodic = true;

            self.xVisualLimits = [0 2 * self.Lx];
            self.yVisualLimits = [-5 * self.L 5 * self.L];
            self.visualScale = 1e3;
            self.name = 'Meandering jet';
        end

        function value = get.k(self)
            % k Return the meander wavenumber.
            value = 2 * pi / self.Lx;
        end

        function thetaValue = theta(self, t, x)
            % theta Evaluate the meander phase.
            thetaValue = self.k * (x - self.cx * t);
        end

        function gammaValue = gamma(self, t, x, y)
            % gamma Evaluate the nondimensional cross-jet coordinate.
            gammaValue = (y - self.A * cos(self.theta(t, x))) / self.L ./ sqrt(1 + (self.k * self.A * sin(self.theta(t, x))).^2);
        end

        function psiValue = psi(self, t, x, y)
            % psi Evaluate the meandering-jet streamfunction.
            psiValue = self.U * self.L * (1 - tanh(self.gamma(t, x, y)));
        end

        function uValue = u(self, t, x, y)
            % u Evaluate the along-jet velocity component.
            uValue = self.U * (1 + (self.k * self.A * sin(self.theta(t, x))).^2).^(-1/2) .* sech(self.gamma(t, x, y)).^2;
        end

        function vValue = v(self, t, x, y)
            % v Evaluate the cross-jet velocity component.
            vValue = -self.U * (self.A * self.k * (1 + self.A * self.A * self.k * self.k) * sin(self.theta(t, x))) .* (1 + (self.k * self.A * sin(self.theta(t, x))).^2).^(-3/2) .* sech(self.gamma(t, x, y)).^2;
        end
    end
end
