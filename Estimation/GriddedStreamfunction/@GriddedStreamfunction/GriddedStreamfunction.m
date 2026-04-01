classdef GriddedStreamfunction < handle
    % Fit a COM-frame streamfunction estimator and observation-event decomposition.
    %
    % `GriddedStreamfunction` fits a mesoscale streamfunction
    % $$\psi(\tilde{x},\tilde{y},t)$$, a center-of-mass trajectory
    % $$m_x(t), m_y(t)$$, and a background velocity
    % $$u^{\mathrm{bg}}(t), v^{\mathrm{bg}}(t)$$ from asynchronous drifter
    % trajectory splines. After the mesoscale fit, the remaining fixed-frame
    % residual is split into background and submesoscale parts at every
    % observation event.
    %
    % The fitted decomposition follows
    %
    % $$
    % \dot{x} = u^{\mathrm{meso}} + u^{\mathrm{bg}} + u^{\mathrm{sm}},
    % \qquad
    % \dot{y} = v^{\mathrm{meso}} + v^{\mathrm{bg}} + v^{\mathrm{sm}},
    % $$
    %
    % with centered coordinates
    % $$\tilde{x} = x - m_x(t)$$ and $$\tilde{y} = y - m_y(t).$$
    %
    % ```matlab
    % fit = GriddedStreamfunction(trajectories);
    % uMeso = fit.uMesoscale(tQuery, xQuery, yQuery);
    % uBg = fit.uBackground(tQuery);
    % ```
    %
    % - Topic: Fit the estimator
    % - Topic: Inspect fitted components
    % - Topic: Inspect decomposition samples
    % - Topic: Evaluate fitted mesoscale
    % - Topic: Evaluate fitted diagnostics
    % - Declaration: classdef GriddedStreamfunction < handle

    properties (SetAccess = private)
        % Fitted COM-frame mesoscale streamfunction spline.
        %
        % This spline evaluates $$\psi(\tilde{x},\tilde{y},t)$$ in the
        % centered coordinates
        % $$\tilde{x} = x - m_x(t)$$ and $$\tilde{y} = y - m_y(t).$$
        %
        % - Topic: Inspect fitted components
        streamfunctionSpline

        % Fitted center-of-mass trajectory.
        %
        % `centerOfMassTrajectory.x(t)` evaluates $$m_x(t)$$ and
        % `centerOfMassTrajectory.y(t)` evaluates $$m_y(t)$$.
        %
        % - Topic: Inspect fitted components
        centerOfMassTrajectory

        % Fitted background-velocity trajectory.
        %
        % `backgroundVelocityTrajectory.x(t)` evaluates
        % $$u^{\mathrm{bg}}(t)$$ and `backgroundVelocityTrajectory.y(t)`
        % evaluates $$v^{\mathrm{bg}}(t)$$.
        %
        % - Topic: Inspect fitted components
        backgroundVelocityTrajectory

        % Fast temporal knot vector used for COM and background fits.
        %
        % - Topic: Inspect fitted components
        fastKnotPoints

        % Fast temporal spline degree for COM and background fits.
        %
        % - Topic: Inspect fitted components
        fastS

        % Mesoscale tensor-product knot vectors `{qKnot, rKnot, tKnot}`.
        %
        % - Topic: Inspect fitted components
        psiKnotPoints

        % Mesoscale spline degrees `[Sq Sr St]`.
        %
        % - Topic: Inspect fitted components
        psiS

        % Representative pooled times from the stride-rule fast basis.
        %
        % This is empty when `fastKnotPoints` are supplied directly.
        %
        % - Topic: Inspect fitted components
        representativeTimes

        % Sorted unique observation times used as trajectory support.
        %
        % - Topic: Inspect fitted components
        fitSupportTimes

        % Number of observations contributed by each drifter.
        %
        % - Topic: Inspect decomposition samples
        sampleCounts

        % Pooled observation times for every drifter event.
        %
        % These times align elementwise with `observedX`, `observedY`,
        % `observedXVelocity`, and `observedYVelocity`.
        %
        % - Topic: Inspect decomposition samples
        observationTimes

        % Pooled observed x-positions.
        %
        % - Topic: Inspect decomposition samples
        observedX

        % Pooled observed y-positions.
        %
        % - Topic: Inspect decomposition samples
        observedY

        % Pooled observed x-velocities.
        %
        % These are the first derivatives of the supplied trajectory
        % splines evaluated at `observationTimes`.
        %
        % - Topic: Inspect decomposition samples
        observedXVelocity

        % Pooled observed y-velocities.
        %
        % These are the first derivatives of the supplied trajectory
        % splines evaluated at `observationTimes`.
        %
        % - Topic: Inspect decomposition samples
        observedYVelocity

        % Pooled COM-frame x-coordinates.
        %
        % - Topic: Inspect decomposition samples
        centeredX

        % Pooled COM-frame y-coordinates.
        %
        % - Topic: Inspect decomposition samples
        centeredY

        % Pooled COM x-velocity samples.
        %
        % - Topic: Inspect decomposition samples
        centerVelocityX

        % Pooled COM y-velocity samples.
        %
        % - Topic: Inspect decomposition samples
        centerVelocityY

        % Pooled mesoscale x-velocity samples.
        %
        % - Topic: Inspect decomposition samples
        uMesoscaleObserved

        % Pooled mesoscale y-velocity samples.
        %
        % - Topic: Inspect decomposition samples
        vMesoscaleObserved

        % Pooled fixed-frame residual x-velocity samples.
        %
        % This is
        % $$\rho^x = \dot{x} - u^{\mathrm{meso}}.$$
        %
        % - Topic: Inspect decomposition samples
        rhoX

        % Pooled fixed-frame residual y-velocity samples.
        %
        % This is
        % $$\rho^y = \dot{y} - v^{\mathrm{meso}}.$$
        %
        % - Topic: Inspect decomposition samples
        rhoY

        % Pooled fitted background x-velocity samples.
        %
        % - Topic: Inspect decomposition samples
        uBackgroundObserved

        % Pooled fitted background y-velocity samples.
        %
        % - Topic: Inspect decomposition samples
        vBackgroundObserved

        % Pooled fitted submesoscale x-velocity samples.
        %
        % - Topic: Inspect decomposition samples
        submesoscaleX

        % Pooled fitted submesoscale y-velocity samples.
        %
        % - Topic: Inspect decomposition samples
        submesoscaleY
    end

    properties (SetAccess = private, Hidden)
        % Developer-only diagnostics from the least-squares fit.
        %
        % This struct stores the design matrices, reduced solves, and
        % fitted coefficient vectors used internally by the estimator.
        %
        % - Topic: Inspect decomposition samples
        % - Developer: true
        fitDiagnostics struct = struct()
    end

    methods
        function self = GriddedStreamfunction(trajectories, options)
            % Fit the estimator from drifter trajectory splines.
            %
            % Use this constructor with one `TrajectorySpline` per drifter.
            % Positions are sampled with `trajectory.x(trajectory.t)` and
            % `trajectory.y(trajectory.t)`, and observed velocities are the
            % first derivatives of the same trajectory splines.
            %
            % The fast temporal basis is used for both the center-of-mass
            % trajectory and the background velocity, while the tensor basis
            % defined by `psiS` and `psiKnotPoints` is used only for the
            % mesoscale streamfunction.
            %
            % - Topic: Fit the estimator
            % - Declaration: self = GriddedStreamfunction(trajectories,psiKnotPoints=...,psiS=...,fastKnotPoints=...,fastS=...)
            % - Parameter trajectories: nonempty vector of `TrajectorySpline` drifters
            % - Parameter psiKnotPoints: optional cell array `{qKnot, rKnot, tKnot}` for the mesoscale basis
            % - Parameter psiS: optional mesoscale spline degree vector `[Sq Sr St]`, default `[2 2 0]`
            % - Parameter fastKnotPoints: optional fast temporal knot vector for COM and background
            % - Parameter fastS: optional fast temporal spline degree, default `3`
            % - Returns self: fitted `GriddedStreamfunction` estimator
            arguments (Input)
                trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
                options.psiKnotPoints = []
                options.psiS (1,3) double {mustBeInteger, mustBeNonnegative} = [2 2 0]
                options.fastKnotPoints = []
                options.fastS (1,1) double {mustBeInteger, mustBeNonnegative} = 3
            end

            self.fitTrajectorySplines(reshape(trajectories, [], 1), ...
                options.psiKnotPoints, reshape(options.psiS, 1, []), ...
                options.fastKnotPoints, options.fastS);
        end

        function values = uBackground(self, t)
            % Evaluate the fitted background x-velocity.
            %
            % - Topic: Evaluate fitted diagnostics
            % - Declaration: values = uBackground(self,t)
            % - Parameter t: evaluation times in seconds
            % - Returns values: background x-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
            end

            values = self.backgroundVelocityTrajectory.x(t);
        end

        function values = vBackground(self, t)
            % Evaluate the fitted background y-velocity.
            %
            % - Topic: Evaluate fitted diagnostics
            % - Declaration: values = vBackground(self,t)
            % - Parameter t: evaluation times in seconds
            % - Returns values: background y-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
            end

            values = self.backgroundVelocityTrajectory.y(t);
        end

        function psiValue = psiMesoscale(self, t, x, y)
            % Evaluate the fitted mesoscale streamfunction.
            %
            % - Topic: Evaluate fitted mesoscale
            % - Declaration: psiValue = psiMesoscale(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns psiValue: mesoscale streamfunction values with the same shape as `x`
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            psiValue = self.streamfunctionSpline.valueAtPoints(q, r, tEval);
        end

        function uValue = uMesoscale(self, t, x, y)
            % Evaluate the mesoscale x-velocity $$-\psi_{\tilde{y}}$$.
            %
            % - Topic: Evaluate fitted mesoscale
            % - Declaration: uValue = uMesoscale(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns uValue: mesoscale x-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            uValue = -self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[0 1 0]);
        end

        function vValue = vMesoscale(self, t, x, y)
            % Evaluate the mesoscale y-velocity $$\psi_{\tilde{x}}$$.
            %
            % - Topic: Evaluate fitted mesoscale
            % - Declaration: vValue = vMesoscale(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns vValue: mesoscale y-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            vValue = self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[1 0 0]);
        end

        function values = sigma_n(self, t, x, y)
            % Evaluate the normal strain field $$\sigma_n = -2\psi_{\tilde{x}\tilde{y}}$$.
            %
            % - Topic: Evaluate fitted diagnostics
            % - Declaration: values = sigma_n(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns values: normal strain in $$s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            values = -2 * self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[1 1 0]);
        end

        function values = sigma_s(self, t, x, y)
            % Evaluate the shear strain field $$\sigma_s = \psi_{\tilde{x}\tilde{x}} - \psi_{\tilde{y}\tilde{y}}$$.
            %
            % - Topic: Evaluate fitted diagnostics
            % - Declaration: values = sigma_s(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns values: shear strain in $$s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            dqq = self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[2 0 0]);
            drr = self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[0 2 0]);
            values = dqq - drr;
        end

        function values = zeta(self, t, x, y)
            % Evaluate the relative-vorticity field $$\zeta = \psi_{\tilde{x}\tilde{x}} + \psi_{\tilde{y}\tilde{y}}$$.
            %
            % - Topic: Evaluate fitted diagnostics
            % - Declaration: values = zeta(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns values: relative vorticity in $$s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, q, r] = centeredCoordinates(self, t, x, y);
            dqq = self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[2 0 0]);
            drr = self.streamfunctionSpline.valueAtPoints(q, r, tEval, D=[0 2 0]);
            values = dqq + drr;
        end
    end

    methods (Access = private)
        fitTrajectorySplines(self, trajectories, psiKnotPoints, psiS, fastKnotPoints, fastS)
        [tEval, q, r] = centeredCoordinates(self, t, x, y)
    end

    methods (Static, Access = private)
        representativeTimes = representativeObservationTimes(tCell)
        fastKnotPoints = validateFastKnotPoints(fastKnotPoints, fastS, allT)
        psiKnotPoints = validatePsiKnotPoints(psiKnotPoints)
        psiKnotPoints = defaultPsiKnotPoints(qAll, rAll, tAll, psiS)
    end
end
