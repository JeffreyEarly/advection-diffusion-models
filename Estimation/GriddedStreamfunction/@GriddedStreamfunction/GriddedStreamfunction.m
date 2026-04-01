classdef GriddedStreamfunction < StreamfunctionModel
    % Fit and evaluate a COM-frame streamfunction on asynchronous drifter trajectories.
    %
    % `GriddedStreamfunction` stores a mesoscale streamfunction
    % $$\psi(\tilde{x},\tilde{y},t)$$ together with trajectory splines for
    % the center of mass
    % $$m_x(t), m_y(t)$$ and the background velocity
    % $$u^{\mathrm{bg}}(t), v^{\mathrm{bg}}(t).$$
    %
    % The total velocity field is
    %
    % $$
    % u(t,x,y) = u^{\mathrm{bg}}(t) - \psi_{\tilde{y}}(\tilde{x},\tilde{y},t),
    % \qquad
    % v(t,x,y) = v^{\mathrm{bg}}(t) + \psi_{\tilde{x}}(\tilde{x},\tilde{y},t),
    % $$
    %
    % with centered coordinates
    % $$\tilde{x} = x - m_x(t)$$ and $$\tilde{y} = y - m_y(t).$$
    %
    % ```matlab
    % fit = GriddedStreamfunction.fitFromTrajectorySplines(trajectories);
    % u = fit.u(tQuery, xQuery, yQuery);
    % ```
    %
    % - Topic: Create the model
    % - Topic: Fit the model
    % - Topic: Inspect fitted components
    % - Topic: Evaluate the streamfunction
    % - Topic: Evaluate the velocity field
    % - Topic: Evaluate fitted diagnostics
    % - Declaration: classdef GriddedStreamfunction < StreamfunctionModel

    properties (SetAccess = private)
        % Centered-frame mesoscale streamfunction spline.
        %
        % This spline is evaluated in the centered coordinates
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
    end

    properties (SetAccess = private, Hidden)
        % Developer-only diagnostics from the fit pipeline.
        %
        % - Topic: Inspect fitted components
        % - Developer: true
        fitState struct = struct()
    end

    methods
        function self = GriddedStreamfunction(streamfunctionSpline, centerOfMassTrajectory, backgroundVelocityTrajectory)
            % Create a gridded streamfunction model from fitted spline components.
            %
            % Use this constructor when the mesoscale streamfunction, the
            % center-of-mass trajectory, and the background velocity have
            % already been fit.
            %
            % - Topic: Create the model
            % - Declaration: self = GriddedStreamfunction(streamfunctionSpline,centerOfMassTrajectory,backgroundVelocityTrajectory)
            % - Parameter streamfunctionSpline: 3-D `TensorSpline` for $$\psi(\tilde{x},\tilde{y},t)$$
            % - Parameter centerOfMassTrajectory: `TrajectorySpline` for $$m_x(t), m_y(t)$$
            % - Parameter backgroundVelocityTrajectory: `TrajectorySpline` for $$u^{\mathrm{bg}}(t), v^{\mathrm{bg}}(t)$$
            % - Returns self: `GriddedStreamfunction` instance
            arguments
                streamfunctionSpline (1,1) TensorSpline
                centerOfMassTrajectory (1,1) TrajectorySpline
                backgroundVelocityTrajectory (1,1) TrajectorySpline
            end

            if streamfunctionSpline.numDimensions ~= 3
                error("GriddedStreamfunction:InvalidStreamfunctionSpline", ...
                    "streamfunctionSpline must be a 3-D TensorSpline with dimensions [q r t].");
            end
            if centerOfMassTrajectory.x.numDimensions ~= 1 || centerOfMassTrajectory.y.numDimensions ~= 1
                error("GriddedStreamfunction:InvalidCenterOfMassTrajectory", ...
                    "centerOfMassTrajectory must store one-dimensional component splines.");
            end
            if backgroundVelocityTrajectory.x.numDimensions ~= 1 || backgroundVelocityTrajectory.y.numDimensions ~= 1
                error("GriddedStreamfunction:InvalidBackgroundVelocityTrajectory", ...
                    "backgroundVelocityTrajectory must store one-dimensional component splines.");
            end

            self.streamfunctionSpline = streamfunctionSpline;
            self.centerOfMassTrajectory = centerOfMassTrajectory;
            self.backgroundVelocityTrajectory = backgroundVelocityTrajectory;
            self.name = "Gridded streamfunction";
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

        function psiValue = psi(self, t, x, y)
            % Evaluate the mesoscale streamfunction on absolute coordinates.
            %
            % - Topic: Evaluate the streamfunction
            % - Declaration: psiValue = psi(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns psiValue: streamfunction values with the same shape as `x`
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
            % - Topic: Evaluate fitted diagnostics
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
            % - Topic: Evaluate fitted diagnostics
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

        function uValue = u(self, t, x, y)
            % Evaluate the total x-velocity.
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: uValue = u(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns uValue: total x-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, ~, ~] = centeredCoordinates(self, t, x, y);
            uValue = self.uBackground(tEval) + self.uMesoscale(tEval, x, y);
        end

        function vValue = v(self, t, x, y)
            % Evaluate the total y-velocity.
            %
            % - Topic: Evaluate the velocity field
            % - Declaration: vValue = v(self,t,x,y)
            % - Parameter t: scalar time or array matching `x` and `y`
            % - Parameter x: x-coordinate array in meters
            % - Parameter y: y-coordinate array in meters
            % - Returns vValue: total y-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
                x {mustBeNumeric, mustBeReal}
                y {mustBeNumeric, mustBeReal}
            end

            [tEval, ~, ~] = centeredCoordinates(self, t, x, y);
            vValue = self.vBackground(tEval) + self.vMesoscale(tEval, x, y);
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

    methods (Static)
        function self = fitFromTrajectorySplines(trajectories, options)
            % Fit the COM, mesoscale streamfunction, and background from drifter trajectory splines.
            %
            % Use this factory with one `TrajectorySpline` per drifter. The
            % fitter samples positions with `trajectory.x(trajectory.t)` and
            % `trajectory.y(trajectory.t)` so it remains compatible with
            % future smoothed trajectory splines.
            %
            % The fast temporal basis is used for both the center-of-mass
            % trajectory and the background velocity, while the slow tensor
            % basis defined by `psiS` and `psiKnotPoints` is used only for
            % the mesoscale streamfunction.
            %
            % - Topic: Fit the model
            % - Declaration: self = fitFromTrajectorySplines(trajectories,options)
            % - Parameter trajectories: nonempty vector of `TrajectorySpline` drifters
            % - Parameter options.psiKnotPoints: optional cell array `{qKnot, rKnot, tKnot}` for the mesoscale basis
            % - Parameter options.psiS: optional mesoscale spline degree vector `[Sq Sr St]`, default `[2 2 0]`
            % - Parameter options.fastKnotPoints: optional fast temporal knot vector for COM and background
            % - Parameter options.fastS: optional fast temporal spline degree, default `3`
            % - Returns self: fitted `GriddedStreamfunction`
            arguments (Input)
                trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
                options.psiKnotPoints = []
                options.psiS (1,3) double {mustBeInteger, mustBeNonnegative} = [2 2 0]
                options.fastKnotPoints = []
                options.fastS (1,1) double {mustBeInteger, mustBeNonnegative} = 3
            end
            arguments (Output)
                self (1,1) GriddedStreamfunction
            end

            trajectories = reshape(trajectories, [], 1);
            nTrajectories = numel(trajectories);
            tCell = cell(nTrajectories, 1);
            xCell = cell(nTrajectories, 1);
            yCell = cell(nTrajectories, 1);
            xDotCell = cell(nTrajectories, 1);
            yDotCell = cell(nTrajectories, 1);

            for iTrajectory = 1:nTrajectories
                ti = trajectories(iTrajectory).t;
                xi = reshape(trajectories(iTrajectory).x(ti), [], 1);
                yi = reshape(trajectories(iTrajectory).y(ti), [], 1);
                xDoti = reshape(trajectories(iTrajectory).x.valueAtPoints(ti, D=1), [], 1);
                yDoti = reshape(trajectories(iTrajectory).y.valueAtPoints(ti, D=1), [], 1);

                tCell{iTrajectory} = ti;
                xCell{iTrajectory} = xi;
                yCell{iTrajectory} = yi;
                xDotCell{iTrajectory} = xDoti;
                yDotCell{iTrajectory} = yDoti;
            end

            psiS = reshape(options.psiS, 1, []);
            allT = vertcat(tCell{:});
            allX = vertcat(xCell{:});
            allY = vertcat(yCell{:});
            allXDot = vertcat(xDotCell{:});
            allYDot = vertcat(yDotCell{:});
            sampleCounts = cellfun(@numel, tCell);
            fitSupportTimes = unique(allT, "sorted");

            if isempty(options.fastKnotPoints)
                representativeTimes = GriddedStreamfunction.representativeObservationTimes(tCell);
                fastSupportTimes = unique(representativeTimes, "sorted");
                if numel(fastSupportTimes) < options.fastS + 1
                    error("GriddedStreamfunction:InsufficientFastSupportTimes", ...
                        "The representative times must contain at least fastS + 1 unique samples.");
                end
                fastKnotPoints = BSpline.knotPointsForDataPoints(fastSupportTimes, S=options.fastS, splineDOF=numel(fastSupportTimes));
            else
                representativeTimes = [];
                fastKnotPoints = GriddedStreamfunction.validateFastKnotPoints(options.fastKnotPoints, options.fastS, allT);
                fastSupportTimes = fitSupportTimes;
            end

            fastBasisMatrix = BSpline.matrixForDataPoints(allT, knotPoints=fastKnotPoints, S=options.fastS);
            fastBasisMatrix = reshape(fastBasisMatrix, numel(allT), []);
            if options.fastS > 0
                fastDerivativeTensor = BSpline.matrixForDataPoints(allT, knotPoints=fastKnotPoints, S=options.fastS, D=1);
                fastDerivativeMatrix = reshape(fastDerivativeTensor(:, :, 2), numel(allT), []);
            else
                fastDerivativeMatrix = zeros(size(fastBasisMatrix));
            end

            centerXCoefficients = fastBasisMatrix \ allX;
            centerYCoefficients = fastBasisMatrix \ allY;
            centerXSpline = TensorSpline(S=options.fastS, knotPoints=fastKnotPoints, xi=centerXCoefficients);
            centerYSpline = TensorSpline(S=options.fastS, knotPoints=fastKnotPoints, xi=centerYCoefficients);
            centerOfMassTrajectory = TrajectorySpline.fromComponentSplines(fitSupportTimes, centerXSpline, centerYSpline);

            mxAll = centerOfMassTrajectory.x(allT);
            myAll = centerOfMassTrajectory.y(allT);
            mxDotAll = fastDerivativeMatrix * centerXCoefficients;
            myDotAll = fastDerivativeMatrix * centerYCoefficients;
            qAll = allX - mxAll;
            rAll = allY - myAll;

            if isempty(options.psiKnotPoints)
                psiKnotPoints = GriddedStreamfunction.defaultPsiKnotPoints(qAll, rAll, allT, psiS);
            else
                psiKnotPoints = GriddedStreamfunction.validatePsiKnotPoints(options.psiKnotPoints);
            end

            try
                psiBasisSize = TensorSpline.basisSizeFromKnotCell(psiKnotPoints, psiS + 1);
            catch
                error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
                    "psiKnotPoints must define at least one basis function in each dimension.");
            end

            qDomain = [psiKnotPoints{1}(1), psiKnotPoints{1}(end)];
            rDomain = [psiKnotPoints{2}(1), psiKnotPoints{2}(end)];
            tDomain = [psiKnotPoints{3}(1), psiKnotPoints{3}(end)];
            if any(qAll < qDomain(1) | qAll > qDomain(2))
                error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
                    "Centered q observations must lie inside the supplied psi spline domain.");
            end
            if any(rAll < rDomain(1) | rAll > rDomain(2))
                error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
                    "Centered r observations must lie inside the supplied psi spline domain.");
            end
            if any(allT < tDomain(1) | allT > tDomain(2))
                error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
                    "Observation times must lie inside the supplied psi time domain.");
            end

            pointMatrix = [qAll, rAll, allT];
            Ru = -sparse(TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[0 1 0]));
            Rv = sparse(TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[1 0 0]));
            designMatrix = [Ru; Rv];
            observedCenteredVelocity = [allXDot - mxDotAll; allYDot - myDotAll];

            numSpatialCoefficients = prod(psiBasisSize(1:2));
            numTimeCoefficients = psiBasisSize(3);
            if numSpatialCoefficients == 1
                reducedSpatialBasis = sparse(1, 0);
                reducedBasisMatrix = sparse(prod(psiBasisSize), 0);
                reducedDesignMatrix = sparse(size(designMatrix, 1), 0);
                reducedCoefficients = zeros(0, 1);
                xi = zeros(prod(psiBasisSize), 1);
            else
                reducedSpatialBasis = spalloc(numSpatialCoefficients, numSpatialCoefficients - 1, ...
                    2 * (numSpatialCoefficients - 1));
                reducedSpatialBasis(1:(numSpatialCoefficients - 1), :) = speye(numSpatialCoefficients - 1);
                reducedSpatialBasis(numSpatialCoefficients, :) = -1;
                reducedBasisMatrix = kron(speye(numTimeCoefficients), reducedSpatialBasis);
                reducedDesignMatrix = designMatrix * reducedBasisMatrix;
                reducedCoefficients = reducedDesignMatrix \ observedCenteredVelocity;
                xi = reducedBasisMatrix * reducedCoefficients;
            end

            streamfunctionSpline = TensorSpline(S=psiS, knotPoints=psiKnotPoints, xi=xi);
            uMesoscaleObserved = -streamfunctionSpline.valueAtPoints(qAll, rAll, allT, D=[0 1 0]);
            vMesoscaleObserved = streamfunctionSpline.valueAtPoints(qAll, rAll, allT, D=[1 0 0]);
            rhoX = allXDot - uMesoscaleObserved;
            rhoY = allYDot - vMesoscaleObserved;

            backgroundXCoefficients = fastBasisMatrix \ rhoX;
            backgroundYCoefficients = fastBasisMatrix \ rhoY;
            backgroundXSpline = TensorSpline(S=options.fastS, knotPoints=fastKnotPoints, xi=backgroundXCoefficients);
            backgroundYSpline = TensorSpline(S=options.fastS, knotPoints=fastKnotPoints, xi=backgroundYCoefficients);
            backgroundVelocityTrajectory = TrajectorySpline.fromComponentSplines(fitSupportTimes, backgroundXSpline, backgroundYSpline);

            uBackgroundObserved = backgroundVelocityTrajectory.x(allT);
            vBackgroundObserved = backgroundVelocityTrajectory.y(allT);
            submesoscaleX = rhoX - uBackgroundObserved;
            submesoscaleY = rhoY - vBackgroundObserved;

            self = GriddedStreamfunction(streamfunctionSpline, centerOfMassTrajectory, backgroundVelocityTrajectory);
            fitState = struct( ...
                "sampleCounts", sampleCounts, ...
                "observationTimes", allT, ...
                "observedX", allX, ...
                "observedY", allY, ...
                "observedXVelocity", allXDot, ...
                "observedYVelocity", allYDot, ...
                "fitSupportTimes", fitSupportTimes, ...
                "representativeTimes", representativeTimes, ...
                "fastS", options.fastS, ...
                "fastKnotPoints", fastKnotPoints, ...
                "fastBasisMatrix", fastBasisMatrix, ...
                "fastDerivativeMatrix", fastDerivativeMatrix, ...
                "centerXCoefficients", centerXCoefficients, ...
                "centerYCoefficients", centerYCoefficients, ...
                "centeredX", qAll, ...
                "centeredY", rAll, ...
                "centerVelocityX", mxDotAll, ...
                "centerVelocityY", myDotAll, ...
                "psiS", psiS, ...
                "psiKnotPoints", {psiKnotPoints}, ...
                "psiBasisSize", psiBasisSize, ...
                "psiReducedSpatialBasis", reducedSpatialBasis, ...
                "psiReducedBasisMatrix", reducedBasisMatrix, ...
                "psiReducedDesignMatrix", reducedDesignMatrix, ...
                "psiDesignMatrix", designMatrix, ...
                "psiObservedCenteredVelocity", observedCenteredVelocity, ...
                "psiReducedCoefficients", reducedCoefficients, ...
                "uMesoscaleObserved", uMesoscaleObserved, ...
                "vMesoscaleObserved", vMesoscaleObserved, ...
                "rhoX", rhoX, ...
                "rhoY", rhoY, ...
                "backgroundXCoefficients", backgroundXCoefficients, ...
                "backgroundYCoefficients", backgroundYCoefficients, ...
                "uBackgroundObserved", uBackgroundObserved, ...
                "vBackgroundObserved", vBackgroundObserved, ...
                "submesoscaleX", submesoscaleX, ...
                "submesoscaleY", submesoscaleY);
            fitState.trajectories = trajectories;
            self.fitState = fitState;
        end
    end

    methods (Access = private)
        function [tEval, q, r] = centeredCoordinates(self, t, x, y)
            if ~isequal(size(x), size(y))
                error("GriddedStreamfunction:EvaluationSizeMismatch", ...
                    "x and y must have the same size.");
            end

            if isscalar(t)
                tEval = repmat(t, size(x));
            elseif isvector(t) && ismatrix(x) && size(x, 1) == numel(t)
                tEval = repmat(reshape(t, [], 1), 1, size(x, 2));
            else
                if ~isequal(size(t), size(x))
                    error("GriddedStreamfunction:EvaluationTimeSizeMismatch", ...
                        "t must be scalar or have the same size as x and y.");
                end
                tEval = t;
            end

            q = x - self.centerOfMassTrajectory.x(tEval);
            r = y - self.centerOfMassTrajectory.y(tEval);
        end
    end

    methods (Static, Access = private)
        function representativeTimes = representativeObservationTimes(tCell)
            pooledTimes = sort(vertcat(tCell{:}));
            representativeTimes = zeros(0, 1);
            pooledIndex = 1;
            while pooledIndex <= numel(pooledTimes)
                tNow = pooledTimes(pooledIndex);
                representativeTimes(end + 1, 1) = tNow; %#ok<AGROW>
                deployedCount = 0;
                for iTrajectory = 1:numel(tCell)
                    ti = tCell{iTrajectory};
                    deployedCount = deployedCount + double(tNow >= ti(1) && tNow <= ti(end));
                end
                if deployedCount < 1
                    error("GriddedStreamfunction:InvalidRepresentativeTimes", ...
                        "The representative-time stride rule produced a zero active-drifter count.");
                end
                pooledIndex = pooledIndex + deployedCount;
            end
        end

        function fastKnotPoints = validateFastKnotPoints(fastKnotPoints, fastS, allT)
            validateattributes(fastKnotPoints, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
            fastKnotPoints = reshape(fastKnotPoints, [], 1);
            if any(diff(fastKnotPoints) < 0)
                error("GriddedStreamfunction:InvalidFastKnotPoints", ...
                    "fastKnotPoints must be nondecreasing.");
            end
            if numel(fastKnotPoints) <= fastS + 1
                error("GriddedStreamfunction:InvalidFastKnotPoints", ...
                    "fastKnotPoints must define at least one fast temporal basis function.");
            end
            if any(allT < fastKnotPoints(1) | allT > fastKnotPoints(end))
                error("GriddedStreamfunction:ObservationOutsideFastDomain", ...
                    "Observation times must lie inside the supplied fast spline domain.");
            end
        end

        function psiKnotPoints = validatePsiKnotPoints(psiKnotPoints)
            if ~iscell(psiKnotPoints) || numel(psiKnotPoints) ~= 3
                error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
                    "psiKnotPoints must be a cell array {qKnot, rKnot, tKnot}.");
            end

            psiKnotPoints = reshape(psiKnotPoints, 1, []);
            for iDim = 1:3
                validateattributes(psiKnotPoints{iDim}, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
                psiKnotPoints{iDim} = reshape(psiKnotPoints{iDim}, [], 1);
                if any(diff(psiKnotPoints{iDim}) < 0)
                    error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
                        "Each psi knot vector must be nondecreasing.");
                end
            end
        end

        function psiKnotPoints = defaultPsiKnotPoints(qAll, rAll, tAll, psiS)
            qMin = min(qAll);
            qMax = max(qAll);
            rMin = min(rAll);
            rMax = max(rAll);
            tMin = min(tAll);
            tMax = max(tAll);

            qPad = max(50, 0.25 * max(qMax - qMin, eps));
            rPad = max(50, 0.25 * max(rMax - rMin, eps));
            tPad = max(1, 0.01 * max(tMax - tMin, eps));

            lowerBounds = [qMin - qPad, rMin - rPad, tMin - tPad];
            upperBounds = [qMax + qPad, rMax + rPad, tMax + tPad];

            psiKnotPoints = cell(1, 3);
            for iDim = 1:3
                repeatCount = psiS(iDim) + 1;
                psiKnotPoints{iDim} = [ ...
                    repmat(lowerBounds(iDim), repeatCount, 1)
                    repmat(upperBounds(iDim), repeatCount, 1)
                    ];
            end
        end
    end
end
