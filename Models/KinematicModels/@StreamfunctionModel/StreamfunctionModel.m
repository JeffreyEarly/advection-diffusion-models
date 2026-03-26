classdef StreamfunctionModel < KinematicModel
    % StreamfunctionModel Base class for incompressible two-dimensional streamfunction models.
    %
    % Streamfunction models define a scalar field $\psi(t,x,y)$ such that
    %
    % $$ u(t,x,y) = -\frac{\partial \psi}{\partial y}, \qquad v(t,x,y) = \frac{\partial \psi}{\partial x}. $$
    %
    % Subclasses implement `psi(self, t, x, y)` and inherit the common
    % streamfunction plotting helper.
    %
    % Topic:
    %   Models

    methods (Abstract)
        % psi Evaluate the streamfunction.
        %
        % Declaration:
        %   `psiValue = psi(self, t, x, y)`
        %
        % Parameters:
        %   `t` - Scalar evaluation time in seconds.
        %
        %   `x`, `y` - Position arrays in meters.
        %
        % Returns:
        %   `psiValue` - Streamfunction values with the same shape as the
        %   input arrays.
        psiValue = psi(self, t, x, y);
    end
end
