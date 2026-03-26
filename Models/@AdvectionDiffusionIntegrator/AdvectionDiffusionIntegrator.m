classdef AdvectionDiffusionIntegrator
    % AdvectionDiffusionIntegrator Integrate particles through a kinematic flow with diffusivity.
    %
    % This class advances particle positions according to
    %
    % $$ d\mathbf{r}(t) = \mathbf{u}(t,\mathbf{r})\,dt + \sqrt{2\kappa}\,d\mathbf{W}_t, $$
    %
    % where $\mathbf{r}(t) = [x(t)\ y(t)]$, $\mathbf{u}(t,\mathbf{r})$ is
    % supplied by a `KinematicModel`, and $\kappa$ is the scalar tracer
    % diffusivity in m^2 s^-1. Boundary limits, periodic wrapping, and
    % polygonal obstacles are delegated to `IntegratorWithObstacles`.
    %
    % The output trajectory arrays follow the contract
    %
    % - `size(t) = [nTimes 1]`
    % - `size(x) = [nTimes nParticles]`
    % - `size(y) = [nTimes nParticles]`
    %
    % Topic:
    %   Models

    properties
        kinematicModel (1,1) KinematicModel
        % Underlying deterministic velocity model.

        kappa (1,1) double {mustBeNonnegative} = 0
        % Scalar diffusivity in m^2 s^-1.

        stepSize (1,1) double {mustBeNonnegative} = 0
        % Internal integrator step size in seconds. A value of `0` uses
        % the requested output increment `dt` as the integration step.
    end

    methods
        function self = AdvectionDiffusionIntegrator(kinematicModel, kappa)
            % AdvectionDiffusionIntegrator Create an advection-diffusion integrator.
            %
            % Declaration:
            %   `self = AdvectionDiffusionIntegrator(kinematicModel, kappa)`
            %
            % Parameters:
            %   `kinematicModel` - `KinematicModel` instance defining
            %   `u(t,x,y)` and `v(t,x,y)`.
            %
            %   `kappa` - Scalar diffusivity in m^2 s^-1.
            arguments
                kinematicModel (1,1) KinematicModel
                kappa (1,1) double {mustBeNonnegative}
            end

            self.kinematicModel = kinematicModel;
            self.kappa = kappa;
        end
    end
end
