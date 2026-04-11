classdef GriddedStreamfunction < handle
    % Fit a COM-frame streamfunction estimator and trajectory decomposition.
    %
    % `GriddedStreamfunction` fits a mesoscale streamfunction
    % $$\psi(\tilde{x},\tilde{y},t)$$, a center-of-mass trajectory
    % $$m_x(t), m_y(t)$$, and an anchored background trajectory
    % $$x^{\mathrm{bg}}(t), y^{\mathrm{bg}}(t)$$ from asynchronous drifter
    % trajectory splines using the rev3 mesoscale-only formulation
    % described in `asynchronous-com-fit-rev3.tex`, with optional hard
    % mesoscale constraints described in `asynchronous-com-fit-rev4.tex`.
    %
    % The estimator first fits the COM position in a fast temporal basis,
    % then solves one least-squares problem for the mesoscale
    % streamfunction coefficients. The background velocity is recovered
    % afterward as the COM residual projected back onto the same fast
    % basis, then integrated to obtain the common anchored background
    % path. When requested, the mesoscale fit imposes a hard
    % `zeroVorticity` or `zeroStrain` constraint on the gauge-reduced
    % streamfunction coefficients.
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
    % The fixed-frame trajectory decomposition is stored as per-drifter
    % spline trajectories satisfying
    %
    % $$
    % x_k = x_k^{\mathrm{bg}} + x_k^{\mathrm{meso}} + x_k^{\mathrm{sm}},
    % \qquad
    % y_k = y_k^{\mathrm{bg}} + y_k^{\mathrm{meso}} + y_k^{\mathrm{sm}},
    % $$
    %
    % where the mesoscale trajectory carries the observed initial drifter
    % position and the background and submesoscale trajectories are
    % zero-anchored at the first drifter sample. In the centered frame,
    % the mesoscale trajectory carries the initial centered position and
    % the centered submesoscale trajectory is zero-anchored.
    %
    % ```matlab
    % fit = GriddedStreamfunction(trajectories);
    % uMeso = fit.uMesoscale(tQuery, xQuery, yQuery);
    % xMeso = fit.decomposition.fixedFrame.mesoscale(1).x(tQuery);
    % decomposition = fit.decomposeTrajectories(otherTrajectories);
    % ```
    %
    % - Topic: Fit the estimator
    % - Topic: Inspect fitted components
    % - Topic: Inspect decomposition trajectories
    % - Topic: Apply fitted decomposition
    % - Topic: Evaluate fitted mesoscale
    % - Topic: Evaluate fitted diagnostics
    % - Topic: Visualize strain angle
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

        % Observed drifter trajectory splines used for the fit.
        %
        % `observedTrajectories` preserves the original
        % `TrajectorySpline` inputs in drifter order.
        %
        % - Topic: Inspect fitted components
        observedTrajectories

        % Fitted center-of-mass trajectory.
        %
        % `centerOfMassTrajectory.x(t)` evaluates $$m_x(t)$$ and
        % `centerOfMassTrajectory.y(t)` evaluates $$m_y(t)$$.
        %
        % - Topic: Inspect fitted components
        centerOfMassTrajectory

        % Fitted common background trajectory.
        %
        % `backgroundTrajectory.x(t)` evaluates $$x^{\mathrm{bg}}(t)$$
        % and `backgroundTrajectory.y(t)` evaluates
        % $$y^{\mathrm{bg}}(t)$$, with
        % $$x^{\mathrm{bg}}(t_0)=0$$ and $$y^{\mathrm{bg}}(t_0)=0$$ at
        % the global fit start time. The recovered background velocity is
        % obtained from `backgroundTrajectory.u(t)` and
        % `backgroundTrajectory.v(t)`.
        %
        % - Topic: Inspect fitted components
        backgroundTrajectory

        % Per-drifter decomposition trajectories in fixed and centered frames.
        %
        % `decomposition.fixedFrame.background`,
        % `decomposition.fixedFrame.mesoscale`, and
        % `decomposition.fixedFrame.submesoscale` are `TrajectorySpline`
        % column vectors aligned with `observedTrajectories`. In the
        % centered frame, `decomposition.centeredFrame.mesoscale` and
        % `decomposition.centeredFrame.submesoscale` store the corresponding
        % COM-frame trajectories.
        %
        % - Topic: Inspect decomposition trajectories
        decomposition struct = struct()

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

        % Hard constraint applied to the fitted mesoscale streamfunction.
        %
        % `mesoscaleConstraint` is `"none"`, `"zeroVorticity"`, or
        % `"zeroStrain"`.
        %
        % - Topic: Inspect fitted components
        mesoscaleConstraint string = "none"

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
            % trajectory and the recovered common background path. The
            % tensor basis defined by `psiS` and `psiKnotPoints` is used
            % only for the mesoscale streamfunction in the centered frame,
            % with only the additive streamfunction gauge removed. Set
            % `mesoscaleConstraint` to impose a hard zero-vorticity or
            % zero-strain mesoscale fit.
            %
            % - Topic: Fit the estimator
            % - Declaration: self = GriddedStreamfunction(trajectories,psiKnotPoints=...,psiS=...,fastKnotPoints=...,fastS=...,mesoscaleConstraint=...)
            % - Parameter trajectories: nonempty vector of `TrajectorySpline` drifters
            % - Parameter psiKnotPoints: optional cell array `{qKnot, rKnot, tKnot}` for the mesoscale basis
            % - Parameter psiS: optional mesoscale spline degree vector `[Sq Sr St]`, default `[2 2 0]`
            % - Parameter fastKnotPoints: optional fast temporal knot vector for COM and background
            % - Parameter fastS: optional fast temporal spline degree, default `3`
            % - Parameter mesoscaleConstraint: optional hard mesoscale constraint `"none"`, `"zeroVorticity"`, or `"zeroStrain"`
            % - Returns self: fitted `GriddedStreamfunction` estimator
            arguments (Input)
                trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
                options.psiKnotPoints = []
                options.psiS (1,3) double {mustBeInteger, mustBeNonnegative} = [2 2 0]
                options.fastKnotPoints = []
                options.fastS (1,1) double {mustBeInteger, mustBeNonnegative} = 3
                options.mesoscaleConstraint {mustBeTextScalar, mustBeMember(options.mesoscaleConstraint, ["none", "zeroVorticity", "zeroStrain"])} = "none"
                options.buildDecomposition (1,1) logical = true
            end

            trajectories = reshape(trajectories, [], 1);
            psiS = reshape(options.psiS, 1, []);
            psiKnotPoints = options.psiKnotPoints;
            fastKnotPoints = options.fastKnotPoints;
            fastS = options.fastS;
            mesoscaleConstraint = string(options.mesoscaleConstraint);

            if ~isempty(fastKnotPoints)
                if ~(isnumeric(fastKnotPoints) && isvector(fastKnotPoints) && ...
                        isreal(fastKnotPoints) && all(isfinite(fastKnotPoints)))
                    error("GriddedStreamfunction:InvalidFastKnotPoints", ...
                        "fastKnotPoints must be a finite real vector.");
                end

                fastKnotPoints = reshape(fastKnotPoints, [], 1);
                if any(diff(fastKnotPoints) < 0)
                    error("GriddedStreamfunction:InvalidFastKnotPoints", ...
                        "fastKnotPoints must be nondecreasing.");
                end
            end

            if ~isempty(psiKnotPoints)
                if ~iscell(psiKnotPoints) || numel(psiKnotPoints) ~= 3
                    error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
                        "psiKnotPoints must be a cell array {qKnot, rKnot, tKnot}.");
                end

                psiKnotPoints = reshape(psiKnotPoints, 1, []);
                for iDim = 1:3
                    knotVector = psiKnotPoints{iDim};
                    if ~(isnumeric(knotVector) && isvector(knotVector) && ...
                            isreal(knotVector) && all(isfinite(knotVector)))
                        error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
                            "Each psi knot vector must be a finite real vector.");
                    end

                    knotVector = reshape(knotVector, [], 1);
                    if any(diff(knotVector) < 0)
                        error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
                            "Each psi knot vector must be nondecreasing.");
                    end

                    psiKnotPoints{iDim} = knotVector;
                end
            end

            fitTrajectorySplines( ...
                self, ...
                trajectories, ...
                psiKnotPoints, ...
                psiS, ...
                fastKnotPoints, ...
                fastS, ...
                mesoscaleConstraint, ...
                options.buildDecomposition);
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

            values = self.backgroundTrajectory.u(t);
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

            values = self.backgroundTrajectory.v(t);
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
        fitTrajectorySplines(self, trajectories, psiKnotPoints, psiS, fastKnotPoints, fastS, mesoscaleConstraint, buildDecomposition)
        [backgroundTrajectory, decomposition] = decomposeTrajectorySet(self, trajectories)
    end

    methods (Access = {?GriddedStreamfunction, ?GriddedStreamfunctionBootstrap, ?GriddedStreamfunctionUnitTests})
        sampleData = decompositionSampleData(self, trajectories)
    end

    methods (Static)
        theta = visualPrincipalStrainAngle(sigma_n, sigma_s, options)
    end

    methods (Static, Access = private)
        representativeTimes = representativeObservationTimes(tCell)
        psiKnotPoints = defaultPsiKnotPoints(qAll, rAll, tAll, psiS)
        Aeq = mesoscaleConstraintMatrix(psiKnotPoints, psiS, G, mesoscaleConstraint)
        fastState = fastBasisState(allT, fastKnotPoints, fastS, shouldComputeDerivative)
    end
end
