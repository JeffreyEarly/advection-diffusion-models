classdef LinearVelocityField < StreamfunctionModel
    % LinearVelocityField Two-dimensional affine velocity field with analytical solutions.
    %
    % This model defines the affine flow
    %
    % $$ \dot{\mathbf{r}} = \mathbf{u}_0 + A\mathbf{r}, $$
    %
    % with
    %
    % $$ \mathbf{u}_0 = \begin{bmatrix} u_0 \\ v_0 \end{bmatrix}, \qquad
    % A = \frac{1}{2}
    % \begin{bmatrix}
    % \sigma_n & \sigma_s - \zeta \\
    % \sigma_s + \zeta & -\sigma_n
    % \end{bmatrix}. $$
    %
    % The strain components satisfy
    %
    % $$ \sigma_n = \sigma \cos(2\theta), \qquad \sigma_s = \sigma \sin(2\theta). $$
    %
    % The method `momentTensorEvolution` integrates the corresponding
    % second-moment system
    %
    % $$ \dot{M} = AM + MA^\top + 2\kappa I. $$
    %
    % Topic:
    %   Models

    properties
        sigma (1,1) double = 0
        % Strain magnitude in s^-1.

        theta (1,1) double = 0
        % Strain orientation in radians.

        zeta (1,1) double = 0
        % Relative vorticity in s^-1.

        u0 (1,1) double = 0
        % Uniform background x-velocity in m s^-1.

        v0 (1,1) double = 0
        % Uniform background y-velocity in m s^-1.
    end

    properties (Dependent)
        sigma_n
        % Normal strain component in s^-1.

        sigma_s
        % Shear strain component in s^-1.
    end

    methods
        function self = LinearVelocityField(sigma, theta, zeta, varargin)
            % LinearVelocityField Create a linear velocity field.
            %
            % Declaration:
            %   `self = LinearVelocityField(sigma, theta, zeta)`
            %   `self = LinearVelocityField(sigma, theta, zeta, u0, v0)`
            %
            % Parameters:
            %   `sigma` - Strain magnitude in s^-1.
            %
            %   `theta` - Strain orientation in radians.
            %
            %   `zeta` - Relative vorticity in s^-1.
            %
            %   `u0`, `v0` - Optional uniform background velocity in
            %   m s^-1. These must be supplied together.
            validateattributes(sigma, {'numeric'}, {'real', 'finite', 'scalar'}, mfilename, 'sigma');
            validateattributes(theta, {'numeric'}, {'real', 'finite', 'scalar'}, mfilename, 'theta');
            validateattributes(zeta, {'numeric'}, {'real', 'finite', 'scalar'}, mfilename, 'zeta');

            if isempty(varargin)
                u0 = 0;
                v0 = 0;
            elseif numel(varargin) == 2
                u0 = varargin{1};
                v0 = varargin{2};
                validateattributes(u0, {'numeric'}, {'real', 'finite', 'scalar'}, mfilename, 'u0');
                validateattributes(v0, {'numeric'}, {'real', 'finite', 'scalar'}, mfilename, 'v0');
            else
                error('LinearVelocityField:InvalidBackgroundVelocity', 'Optional arguments `u0` and `v0` must be supplied together.');
            end

            self.xVisualLimits = 1e3 * [-1 1];
            self.yVisualLimits = 1e3 * [-1 1];
            self.u0 = u0;
            self.v0 = v0;
            self.zeta = zeta;
            self.sigma = sigma;
            self.theta = theta;
            self.updateName();
        end

        function set.zeta(self, value)
            self.zeta = value;
            self.updateName();
        end

        function set.sigma(self, value)
            self.sigma = value;
            self.updateName();
        end

        function set.theta(self, value)
            self.theta = value;
            self.updateName();
        end

        function value = get.sigma_n(self)
            value = self.sigma * cos(2 * self.theta);
        end

        function value = get.sigma_s(self)
            value = self.sigma * sin(2 * self.theta);
        end

        function psiValue = psi(self, t, x, y)
            % psi Evaluate the quadratic streamfunction.
            %
            % Declaration:
            %   `psiValue = psi(self, t, x, y)`
            %
            % The implemented streamfunction is
            %
            % $$ \psi = -u_0 y + v_0 x + \frac{1}{4}(\sigma_s + \zeta)x^2 - \frac{1}{2}\sigma_n xy - \frac{1}{4}(\sigma_s - \zeta)y^2. $$
            psiValue = -self.u0 * y + self.v0 * x + 0.25 * (self.sigma_s + self.zeta) * x .* x - 0.5 * self.sigma_n * x .* y - 0.25 * (self.sigma_s - self.zeta) * y .* y;
        end

        function uValue = u(self, t, x, y)
            % u Evaluate the affine x-velocity.
            %
            % Declaration:
            %   `uValue = u(self, t, x, y)`
            %
            % The implemented velocity component is
            %
            % $$ u = u_0 + \frac{1}{2}\left(\sigma_n x + (\sigma_s - \zeta)y\right). $$
            uValue = self.u0 + 0.5 * (self.sigma_n * x + (self.sigma_s - self.zeta) * y);
        end

        function vValue = v(self, t, x, y)
            % v Evaluate the affine y-velocity.
            %
            % Declaration:
            %   `vValue = v(self, t, x, y)`
            %
            % The implemented velocity component is
            %
            % $$ v = v_0 + \frac{1}{2}\left((\sigma_s + \zeta)x - \sigma_n y\right). $$
            vValue = self.v0 + 0.5 * ((self.sigma_s + self.zeta) * x - self.sigma_n * y);
        end
    end

    methods (Static)
        function [sigma_n, sigma_s] = normalAndShearFromSigmaTheta(sigma, theta)
            % normalAndShearFromSigmaTheta Convert strain magnitude and orientation to components.
            %
            % Declaration:
            %   `[sigma_n, sigma_s] = normalAndShearFromSigmaTheta(sigma, theta)`
            sigma_n = sigma * cos(2 * theta);
            sigma_s = sigma * sin(2 * theta);
        end

        function [sigma, theta] = sigmaThetaFromNormalAndShear(sigma_n, sigma_s)
            % sigmaThetaFromNormalAndShear Convert normal and shear strain to magnitude and orientation.
            %
            % Declaration:
            %   `[sigma, theta] = sigmaThetaFromNormalAndShear(sigma_n, sigma_s)`
            sigma = sqrt(sigma_n .* sigma_n + sigma_s .* sigma_s);
            theta = atan2(sigma_s, sigma_n) / 2;
        end
    end

    methods (Access = private)
        function updateName(self)
            if self.sigma * self.sigma > self.zeta * self.zeta
                self.name = 'Strain-dominated linear velocity field';
            elseif self.sigma * self.sigma < self.zeta * self.zeta
                self.name = 'Vorticity-dominated linear velocity field';
            else
                self.name = 'Strain-vorticity matched linear velocity field';
            end
        end
    end
end
