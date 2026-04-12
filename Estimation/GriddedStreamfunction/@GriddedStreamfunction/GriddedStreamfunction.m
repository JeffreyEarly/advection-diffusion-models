classdef GriddedStreamfunction < CAAnnotatedClass
    % Fit a COM-frame streamfunction estimator and trajectory decomposition.
    %
    % `GriddedStreamfunction` fits a mesoscale streamfunction
    % $$\psi(\tilde{x},\tilde{y},t)$$, a center-of-mass trajectory
    % $$m_x(t), m_y(t)$$, and an anchored background trajectory
    % $$x^{\mathrm{bg}}(t), y^{\mathrm{bg}}(t)$$ from asynchronous drifter
    % trajectory splines using the coupled estimator described in
    % `asynchronous-com-fit-rev6.tex`.
    %
    % The estimator uses the fast temporal spline basis both to define
    % the COM smoothing operator and to represent the common background
    % velocity. It first fits the COM trajectory, then solves one
    % least-squares problem for the gauge-reduced mesoscale
    % streamfunction coefficients in the COM frame, and finally recovers
    % the background velocity from the COM-space residual in that same
    % fast basis. Optional hard `zeroVorticity` and `zeroStrain`
    % constraints are applied directly to the mesoscale coefficients.
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
    % fit = GriddedStreamfunction.fromTrajectories(trajectories);
    % uMeso = fit.uMesoscale(tQuery, xQuery, yQuery);
    % xMeso = fit.decomposition.fixedFrame.mesoscale(1).x(tQuery);
    % decomposition = fit.decomposeTrajectories(otherTrajectories);
    % ```
    %
    % - Topic: Fit the estimator
    % - Topic: Read from file
    % - Topic: Write to file
    % - Topic: Inspect fitted components
    % - Topic: Inspect decomposition trajectories
    % - Topic: Apply fitted decomposition
    % - Topic: Evaluate fitted mesoscale
    % - Topic: Evaluate fitted diagnostics
    % - Topic: Visualize strain angle
    % - Declaration: classdef GriddedStreamfunction < CAAnnotatedClass

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
        % `backgroundTrajectory.v(t)`. This is the single shared
        % background path for the fitted estimator; for drifter `k`,
        % `decomposition.fixedFrame.background(k)` is the same path
        % re-anchored so it starts at zero at that drifter's first sample
        % time.
        %
        % - Topic: Inspect fitted components
        backgroundTrajectory

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

        % Identifiable mesoscale degrees of freedom after gauge and constraints.
        %
        % `mesoscaleDegreesOfFreedom` counts the solved mesoscale spline
        % degrees of freedom after removing the additive streamfunction
        % gauge and applying the selected hard mesoscale constraint. This
        % scalar is determined by the resolved mesoscale basis and the
        % chosen structural constraint, not by the fitted coefficient
        % values themselves.
        %
        % - Topic: Inspect fitted components
        mesoscaleDegreesOfFreedom (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    end

    properties (Dependent)
        % Per-drifter decomposition trajectories in fixed and centered frames.
        %
        % Use `fit.decomposition` to inspect how the fitted estimator
        % reconstructs the same drifters that were used in the fit. Each
        % field is a `TrajectorySpline` column vector aligned one-for-one
        % with `observedTrajectories`, so
        % `decomposition.fixedFrame.mesoscale(k)` is the mesoscale
        % trajectory for `observedTrajectories(k)`. Call
        % `decomposeTrajectories` when the same fitted COM and mesoscale
        % fields should be applied to a different set of drifters.
        %
        % In the fixed frame,
        % `decomposition.fixedFrame.background`,
        % `decomposition.fixedFrame.mesoscale`, and
        % `decomposition.fixedFrame.submesoscale` satisfy
        %
        % $$
        % x_k = x_k^{\mathrm{bg}} + x_k^{\mathrm{meso}} + x_k^{\mathrm{sm}},
        % \qquad
        % y_k = y_k^{\mathrm{bg}} + y_k^{\mathrm{meso}} + y_k^{\mathrm{sm}}.
        % $$
        %
        % The canonical shared background path itself is stored in
        % `backgroundTrajectory`. The fixed-frame background component is
        % that same path re-anchored so that
        % $$x_k^{\mathrm{bg}}(t_{k,0}) = y_k^{\mathrm{bg}}(t_{k,0}) = 0$$
        % at the first sample time $$t_{k,0}$$ of drifter $$k$$. The
        % fixed-frame mesoscale trajectory carries the observed initial
        % drifter position, and the fixed-frame submesoscale residual is
        % zero-anchored.
        %
        % In the centered frame,
        % `decomposition.centeredFrame.mesoscale` and
        % `decomposition.centeredFrame.submesoscale` satisfy
        %
        % $$
        % \tilde{x}_k = \tilde{x}_k^{\mathrm{meso}} + \tilde{x}_k^{\mathrm{sm}},
        % \qquad
        % \tilde{y}_k = \tilde{y}_k^{\mathrm{meso}} + \tilde{y}_k^{\mathrm{sm}},
        % $$
        %
        % after subtracting `centerOfMassTrajectory` from the observed
        % drifter path. Use these centered-frame splines when you want to
        % compare the fitted mesoscale motion against the observed
        % COM-relative motion.
        %
        % ```matlab
        % iDrifter = 1;
        % trajectory = fit.observedTrajectories(iDrifter);
        % ti = trajectory.t;
        %
        % background = fit.decomposition.fixedFrame.background(iDrifter);
        % mesoscale = fit.decomposition.fixedFrame.mesoscale(iDrifter);
        % submesoscale = fit.decomposition.fixedFrame.submesoscale(iDrifter);
        % xRecon = background.x(ti) + mesoscale.x(ti) + submesoscale.x(ti);
        %
        % centeredMesoscale = fit.decomposition.centeredFrame.mesoscale(iDrifter);
        % centeredSubmesoscale = fit.decomposition.centeredFrame.submesoscale(iDrifter);
        % [~, qObs, rObs] = fit.centeredCoordinates(ti, trajectory.x(ti), trajectory.y(ti));
        % qRecon = centeredMesoscale.x(ti) + centeredSubmesoscale.x(ti);
        % rRecon = centeredMesoscale.y(ti) + centeredSubmesoscale.y(ti);
        % ```
        %
        % - Topic: Inspect decomposition trajectories
        decomposition

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
    end

    properties (Hidden, SetAccess = private)
        fixedFrameBackgroundTrajectories = TrajectorySpline.empty(0, 1)
        fixedFrameMesoscaleTrajectories = TrajectorySpline.empty(0, 1)
        fixedFrameSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1)
        centeredFrameMesoscaleTrajectories = TrajectorySpline.empty(0, 1)
        centeredFrameSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1)
    end

    properties (Access = private)
        observedTrajectorySampleData struct = struct([])
    end

    methods
        function self = GriddedStreamfunction(options)
            % Create a fit from canonical solved-state properties.
            %
            % Use this low-level constructor when the fitted spline state
            % is already available, for example after reading a persisted
            % restart file. For fitting from observed drifter trajectories,
            % use `GriddedStreamfunction.fromTrajectories(...)`.
            %
            % ```matlab
            % fit = GriddedStreamfunction( ...
            %     streamfunctionSpline=psiSpline, ...
            %     observedTrajectories=trajectories, ...
            %     centerOfMassTrajectory=centerTrajectory, ...
            %     backgroundTrajectory=backgroundTrajectory, ...
            %     fixedFrameBackgroundTrajectories=backgroundTrajectories, ...
            %     fixedFrameMesoscaleTrajectories=mesoscaleTrajectories, ...
            %     fixedFrameSubmesoscaleTrajectories=submesoscaleTrajectories, ...
            %     centeredFrameMesoscaleTrajectories=centeredMesoscaleTrajectories, ...
            %     centeredFrameSubmesoscaleTrajectories=centeredSubmesoscaleTrajectories, ...
            %     mesoscaleConstraint="none", ...
            %     representativeTimes=representativeTimes, ...
            %     fitSupportTimes=fitSupportTimes);
            % ```
            %
            % - Topic: Fit the estimator
            % - Declaration: self = GriddedStreamfunction(options)
            % - Parameter options.streamfunctionSpline: fitted centered-frame mesoscale streamfunction spline
            % - Parameter options.observedTrajectories: observed drifter trajectory splines aligned with the fit
            % - Parameter options.centerOfMassTrajectory: fitted center-of-mass trajectory
            % - Parameter options.backgroundTrajectory: fitted common background trajectory
            % - Parameter options.fixedFrameBackgroundTrajectories: stored fixed-frame background decomposition trajectories
            % - Parameter options.fixedFrameMesoscaleTrajectories: stored fixed-frame mesoscale decomposition trajectories
            % - Parameter options.fixedFrameSubmesoscaleTrajectories: stored fixed-frame submesoscale decomposition trajectories
            % - Parameter options.centeredFrameMesoscaleTrajectories: stored centered-frame mesoscale decomposition trajectories
            % - Parameter options.centeredFrameSubmesoscaleTrajectories: stored centered-frame submesoscale decomposition trajectories
            % - Parameter options.mesoscaleConstraint: hard mesoscale constraint `"none"`, `"zeroVorticity"`, or `"zeroStrain"`
            % - Parameter options.representativeTimes: pooled representative times from the stride-rule fast basis
            % - Parameter options.fitSupportTimes: sorted unique observation times used as trajectory support
            % - Parameter options.mesoscaleDegreesOfFreedom: stored identifiable mesoscale degrees of freedom after gauge reduction and hard constraints
            % - Returns self: canonical `GriddedStreamfunction` instance
            arguments
                options.streamfunctionSpline = TensorSpline.empty(0, 1)
                options.observedTrajectories = TrajectorySpline.empty(0, 1)
                options.centerOfMassTrajectory = TrajectorySpline.empty(0, 1)
                options.backgroundTrajectory = TrajectorySpline.empty(0, 1)
                options.fixedFrameBackgroundTrajectories = TrajectorySpline.empty(0, 1)
                options.fixedFrameMesoscaleTrajectories = TrajectorySpline.empty(0, 1)
                options.fixedFrameSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1)
                options.centeredFrameMesoscaleTrajectories = TrajectorySpline.empty(0, 1)
                options.centeredFrameSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1)
                options.mesoscaleConstraint {mustBeTextScalar, mustBeMember(options.mesoscaleConstraint, ["none", "zeroVorticity", "zeroStrain"])} = "none"
                options.representativeTimes (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.fitSupportTimes (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.mesoscaleDegreesOfFreedom (1,1) double {mustBeInteger,mustBeNonnegative} = 0
            end

            self@CAAnnotatedClass();

            if nargin == 0
                return
            end

            GriddedStreamfunction.validateCanonicalState(options);
            self.streamfunctionSpline = options.streamfunctionSpline;
            self.observedTrajectories = reshape(options.observedTrajectories, [], 1);
            self.centerOfMassTrajectory = options.centerOfMassTrajectory;
            self.backgroundTrajectory = reshape(options.backgroundTrajectory, [], 1);
            self.fixedFrameBackgroundTrajectories = reshape(options.fixedFrameBackgroundTrajectories, [], 1);
            self.fixedFrameMesoscaleTrajectories = reshape(options.fixedFrameMesoscaleTrajectories, [], 1);
            self.fixedFrameSubmesoscaleTrajectories = reshape(options.fixedFrameSubmesoscaleTrajectories, [], 1);
            self.centeredFrameMesoscaleTrajectories = reshape(options.centeredFrameMesoscaleTrajectories, [], 1);
            self.centeredFrameSubmesoscaleTrajectories = reshape(options.centeredFrameSubmesoscaleTrajectories, [], 1);
            self.mesoscaleConstraint = string(options.mesoscaleConstraint);
            self.representativeTimes = options.representativeTimes;
            self.fitSupportTimes = options.fitSupportTimes;
            self.mesoscaleDegreesOfFreedom = options.mesoscaleDegreesOfFreedom;
            self.refreshObservedTrajectorySampleData();
        end
    end

    methods (Static)
        function self = fromFile(path)
            % Read a fitted estimator from a NetCDF restart file.
            %
            % `fromFile` reconstructs the canonical solved state written by
            % `writeToFile` without rerunning the trajectory fit.
            %
            % - Topic: Read from file
            % - Declaration: self = fromFile(path)
            % - Parameter path: path to the NetCDF restart file
            % - Returns self: reconstructed `GriddedStreamfunction` estimator
            arguments (Input)
                path {mustBeTextScalar}
            end

            self = griddedStreamfunctionFromFile(char(path));
        end

        function self = fromTrajectories(trajectories, options)
            % Fit the estimator from drifter trajectory splines.
            %
            % Use this factory with one `TrajectorySpline` per drifter.
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
            % - Declaration: self = fromTrajectories(trajectories,psiKnotPoints=...,psiS=...,fastKnotPoints=...,fastS=...,mesoscaleConstraint=...)
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
                options.sampleData struct = struct([])
            end

            trajectories = reshape(trajectories, [], 1);
            psiS = reshape(options.psiS, 1, []);
            psiKnotPoints = options.psiKnotPoints;
            fastKnotPoints = options.fastKnotPoints;
            fastS = options.fastS;
            mesoscaleConstraint = string(options.mesoscaleConstraint);
            sampleData = options.sampleData;

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

            self = GriddedStreamfunction();
            fitTrajectorySplines( ...
                self, ...
                trajectories, ...
                psiKnotPoints, ...
                psiS, ...
                fastKnotPoints, ...
                fastS, ...
                mesoscaleConstraint, ...
                options.buildDecomposition, ...
                sampleData);
        end
    end

    methods
        function decomposition = get.decomposition(self)
            if isempty(self.fixedFrameBackgroundTrajectories)
                decomposition = struct();
                return
            end

            decomposition = struct( ...
                "fixedFrame", struct( ...
                    "background", self.fixedFrameBackgroundTrajectories, ...
                    "mesoscale", self.fixedFrameMesoscaleTrajectories, ...
                    "submesoscale", self.fixedFrameSubmesoscaleTrajectories), ...
                "centeredFrame", struct( ...
                    "mesoscale", self.centeredFrameMesoscaleTrajectories, ...
                    "submesoscale", self.centeredFrameSubmesoscaleTrajectories));
        end

        function fastKnotPoints = get.fastKnotPoints(self)
            if isempty(self.centerOfMassTrajectory)
                fastKnotPoints = zeros(0, 1);
                return
            end

            fastKnotPoints = reshape(self.centerOfMassTrajectory.x.knotPoints, [], 1);
        end

        function fastS = get.fastS(self)
            if isempty(self.centerOfMassTrajectory)
                fastS = [];
                return
            end

            fastS = reshape(self.centerOfMassTrajectory.x.S, 1, []);
        end

        function psiKnotPoints = get.psiKnotPoints(self)
            if isempty(self.streamfunctionSpline)
                psiKnotPoints = cell(1, 0);
                return
            end

            psiKnotPoints = self.streamfunctionSpline.knotPoints;
        end

        function psiS = get.psiS(self)
            if isempty(self.streamfunctionSpline)
                psiS = [];
                return
            end

            psiS = reshape(self.streamfunctionSpline.S, 1, []);
        end
    end

    methods (Access = protected)
        function restoreOptionalPersistedPropertiesFromGroup(self, ~)
            self.refreshObservedTrajectorySampleData();
        end
    end

    methods
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
        fitTrajectorySplines(self, trajectories, psiKnotPoints, psiS, fastKnotPoints, fastS, mesoscaleConstraint, buildDecomposition, sampleData)
        [backgroundTrajectory, decomposition] = decomposeTrajectorySet(self, trajectories)
    end

    methods (Access = {?GriddedStreamfunction, ?GriddedStreamfunctionBootstrap, ?GriddedStreamfunctionUnitTests})
        sampleData = decompositionSampleData(self, trajectories)
    end

    methods (Static)
        theta = visualPrincipalStrainAngle(sigma_n, sigma_s, options)
    end

    methods (Static, Access = {?GriddedStreamfunction, ?GriddedStreamfunctionBootstrap, ?GriddedStreamfunctionUnitTests, ?GriddedStreamfunctionBootstrapUnitTests})
        sampleData = sampleTrajectoryData(trajectories)
        sampleData = resampledTrajectorySampleData(observedSampleData, sampledIndices)
    end

    methods (Static, Hidden)
        function self = annotatedClassFromGroup(group)
            vars = CAAnnotatedClass.propertyValuesFromGroup(group, { ...
                'observedTrajectories', ...
                'fixedFrameBackgroundTrajectories', ...
                'fixedFrameMesoscaleTrajectories', ...
                'fixedFrameSubmesoscaleTrajectories', ...
                'centeredFrameMesoscaleTrajectories', ...
                'centeredFrameSubmesoscaleTrajectories', ...
                'mesoscaleConstraint', ...
                'representativeTimes', ...
                'fitSupportTimes', ...
                'mesoscaleDegreesOfFreedom'});
            vars.streamfunctionSpline = TensorSpline.annotatedClassFromGroup(group.groupWithName('streamfunctionSpline'));
            vars.centerOfMassTrajectory = TrajectorySpline.annotatedClassFromGroup(group.groupWithName('centerOfMassTrajectory'));
            vars.backgroundTrajectory = TrajectorySpline.annotatedClassFromGroup(group.groupWithName('backgroundTrajectory'));
            GriddedStreamfunction.validateCanonicalState(vars);

            self = GriddedStreamfunction();
            self.streamfunctionSpline = vars.streamfunctionSpline;
            self.observedTrajectories = reshape(vars.observedTrajectories, [], 1);
            self.centerOfMassTrajectory = vars.centerOfMassTrajectory;
            self.backgroundTrajectory = reshape(vars.backgroundTrajectory, [], 1);
            self.fixedFrameBackgroundTrajectories = reshape(vars.fixedFrameBackgroundTrajectories, [], 1);
            self.fixedFrameMesoscaleTrajectories = reshape(vars.fixedFrameMesoscaleTrajectories, [], 1);
            self.fixedFrameSubmesoscaleTrajectories = reshape(vars.fixedFrameSubmesoscaleTrajectories, [], 1);
            self.centeredFrameMesoscaleTrajectories = reshape(vars.centeredFrameMesoscaleTrajectories, [], 1);
            self.centeredFrameSubmesoscaleTrajectories = reshape(vars.centeredFrameSubmesoscaleTrajectories, [], 1);
            self.mesoscaleConstraint = string(vars.mesoscaleConstraint);
            self.representativeTimes = vars.representativeTimes;
            self.fitSupportTimes = vars.fitSupportTimes;
            self.mesoscaleDegreesOfFreedom = vars.mesoscaleDegreesOfFreedom;
            self.refreshObservedTrajectorySampleData();
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0, 0);
            propertyAnnotations(end+1) = CAObjectProperty('streamfunctionSpline', 'Fitted centered-frame mesoscale streamfunction spline.');
            propertyAnnotations(end+1) = CAObjectProperty('observedTrajectories', 'Observed drifter trajectory splines used for the fit.');
            propertyAnnotations(end+1) = CAObjectProperty('centerOfMassTrajectory', 'Fitted center-of-mass trajectory.');
            propertyAnnotations(end+1) = CAObjectProperty('backgroundTrajectory', 'Fitted common background trajectory.');
            propertyAnnotations(end+1) = CAObjectProperty('fixedFrameBackgroundTrajectories', 'Stored fixed-frame background decomposition trajectories.');
            propertyAnnotations(end+1) = CAObjectProperty('fixedFrameMesoscaleTrajectories', 'Stored fixed-frame mesoscale decomposition trajectories.');
            propertyAnnotations(end+1) = CAObjectProperty('fixedFrameSubmesoscaleTrajectories', 'Stored fixed-frame submesoscale decomposition trajectories.');
            propertyAnnotations(end+1) = CAObjectProperty('centeredFrameMesoscaleTrajectories', 'Stored centered-frame mesoscale decomposition trajectories.');
            propertyAnnotations(end+1) = CAObjectProperty('centeredFrameSubmesoscaleTrajectories', 'Stored centered-frame submesoscale decomposition trajectories.');
            propertyAnnotations(end+1) = CAPropertyAnnotation('mesoscaleConstraint', 'Hard mesoscale constraint applied to the fit.');
            propertyAnnotations(end+1) = CADimensionProperty('representativeTimes', '', 'Representative pooled times from the stride-rule fast basis.');
            propertyAnnotations(end+1) = CADimensionProperty('fitSupportTimes', '', 'Sorted unique observation times used as trajectory support.');
            propertyAnnotations(end+1) = CANumericProperty('mesoscaleDegreesOfFreedom', {}, '', 'Stored identifiable mesoscale degrees of freedom after gauge reduction and hard constraints.');
        end

        function names = classRequiredPropertyNames()
            names = { ...
                'streamfunctionSpline', ...
                'observedTrajectories', ...
                'centerOfMassTrajectory', ...
                'backgroundTrajectory', ...
                'fixedFrameBackgroundTrajectories', ...
                'fixedFrameMesoscaleTrajectories', ...
                'fixedFrameSubmesoscaleTrajectories', ...
                'centeredFrameMesoscaleTrajectories', ...
                'centeredFrameSubmesoscaleTrajectories', ...
                'mesoscaleConstraint', ...
                'representativeTimes', ...
                'fitSupportTimes', ...
                'mesoscaleDegreesOfFreedom'};
        end
    end

    methods (Static, Access = private)
        function validateCanonicalState(options)
            if isempty(options.streamfunctionSpline) && isempty(options.observedTrajectories) && ...
                    isempty(options.centerOfMassTrajectory) && isempty(options.backgroundTrajectory) && ...
                    isempty(options.fixedFrameBackgroundTrajectories) && isempty(options.fixedFrameMesoscaleTrajectories) && ...
                    isempty(options.fixedFrameSubmesoscaleTrajectories) && isempty(options.centeredFrameMesoscaleTrajectories) && ...
                    isempty(options.centeredFrameSubmesoscaleTrajectories) && isempty(options.representativeTimes) && ...
                    isempty(options.fitSupportTimes) && string(options.mesoscaleConstraint) == "none" && ...
                    options.mesoscaleDegreesOfFreedom == 0
                return
            end

            if ~(isa(options.streamfunctionSpline, 'TensorSpline') && isscalar(options.streamfunctionSpline))
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'streamfunctionSpline must be a scalar TensorSpline.');
            end
            if ~(isa(options.observedTrajectories, 'TrajectorySpline') && isvector(options.observedTrajectories) && ~isempty(options.observedTrajectories))
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'observedTrajectories must be a nonempty vector of TrajectorySpline objects.');
            end
            if ~(isa(options.centerOfMassTrajectory, 'TrajectorySpline') && isscalar(options.centerOfMassTrajectory))
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'centerOfMassTrajectory must be a scalar TrajectorySpline.');
            end
            if ~(isa(options.backgroundTrajectory, 'TrajectorySpline') && isvector(options.backgroundTrajectory))
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'backgroundTrajectory must be a TrajectorySpline or empty TrajectorySpline array.');
            end

            trajectoryArrays = { ...
                options.fixedFrameBackgroundTrajectories, ...
                options.fixedFrameMesoscaleTrajectories, ...
                options.fixedFrameSubmesoscaleTrajectories, ...
                options.centeredFrameMesoscaleTrajectories, ...
                options.centeredFrameSubmesoscaleTrajectories};
            for iArray = 1:numel(trajectoryArrays)
                if ~(isa(trajectoryArrays{iArray}, 'TrajectorySpline') && isvector(trajectoryArrays{iArray}))
                    error('GriddedStreamfunction:InvalidCanonicalState', ...
                        'Stored decomposition trajectories must be TrajectorySpline vectors.');
                end
            end

            if ~isempty(options.backgroundTrajectory) && ~isscalar(options.backgroundTrajectory)
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'backgroundTrajectory must be a scalar TrajectorySpline when it is present.');
            end
            if ~isscalar(options.mesoscaleDegreesOfFreedom) || ~isfinite(options.mesoscaleDegreesOfFreedom) || ...
                    options.mesoscaleDegreesOfFreedom < 0 || round(options.mesoscaleDegreesOfFreedom) ~= options.mesoscaleDegreesOfFreedom
                error('GriddedStreamfunction:InvalidCanonicalState', ...
                    'mesoscaleDegreesOfFreedom must be a finite nonnegative integer scalar.');
            end
        end

        representativeTimes = representativeObservationTimes(tCell)
        psiKnotPoints = defaultPsiKnotPoints(qAll, rAll, tAll, psiS)
        Aeq = mesoscaleConstraintMatrix(psiKnotPoints, psiS, G, mesoscaleConstraint)
        fastState = fastBasisState(allT, fastKnotPoints, fastS, shouldComputeDerivative)
    end

    methods (Access = private)
        function refreshObservedTrajectorySampleData(self)
            if isempty(self.observedTrajectories)
                self.observedTrajectorySampleData = struct([]);
            else
                self.observedTrajectorySampleData = GriddedStreamfunction.sampleTrajectoryData(self.observedTrajectories);
            end
        end
    end
end
