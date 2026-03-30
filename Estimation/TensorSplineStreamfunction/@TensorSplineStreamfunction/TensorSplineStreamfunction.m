classdef TensorSplineStreamfunction < StreamfunctionModel
    % Tensor-spline streamfunction fit in a moving center-of-mass frame.
    %
    % `TensorSplineStreamfunction` stores a streamfunction spline in the
    % centered coordinates
    %
    % $$
    % q = x - x_c(t), \qquad r = y - y_c(t),
    % $$
    %
    % together with 1-D splines for the center trajectory $$x_c(t)$$ and
    % $$y_c(t)$$. The total velocity field is then
    %
    % $$
    % u(t,x,y) = \frac{d x_c}{dt} - \frac{\partial \psi}{\partial r}(q,r,t),
    % \qquad
    % v(t,x,y) = \frac{d y_c}{dt} + \frac{\partial \psi}{\partial q}(q,r,t).
    % $$
    %
    % Use `fitFromTrajectories(...)` for asynchronous drifter cell arrays
    % or `fitFromSynchronousTrajectories(...)` when all drifters share a
    % common observation-time vector.
    %
    % ```matlab
    % model = TensorSplineStreamfunction.fitFromTrajectories(xCell, yCell, tCell);
    % u = model.u(tQuery, xQuery, yQuery);
    % ```
    %
    % - Topic: Create the model
    % - Topic: Fit the model
    % - Topic: Inspect fitted components
    % - Topic: Evaluate the streamfunction
    % - Topic: Evaluate the velocity field
    % - Topic: Evaluate fitted diagnostics
    % - Declaration: classdef TensorSplineStreamfunction < StreamfunctionModel

    properties (SetAccess = private)
        % Centered-frame streamfunction spline $$\psi(q,r,t)$$.
        %
        % This 3-D spline is evaluated on the centered coordinates
        % `q = x - xCenter(t)` and `r = y - yCenter(t)`.
        %
        % - Topic: Inspect fitted components
        streamfunctionSpline

        % Center-trajectory spline $$x_c(t)$$ in meters.
        %
        % - Topic: Inspect fitted components
        centerXSpline

        % Center-trajectory spline $$y_c(t)$$ in meters.
        %
        % - Topic: Inspect fitted components
        centerYSpline
    end

    properties (SetAccess = private, Hidden)
        fitState struct = struct()
    end

    methods
        function self = TensorSplineStreamfunction(streamfunctionSpline, centerXSpline, centerYSpline)
            % Create a fitted tensor-spline streamfunction model.
            %
            % Pass a 3-D streamfunction spline in centered coordinates and
            % 1-D splines for the center trajectory components.
            %
            % - Topic: Create the model
            % - Declaration: self = TensorSplineStreamfunction(streamfunctionSpline,centerXSpline,centerYSpline)
            % - Parameter streamfunctionSpline: 3-D `TensorSpline` for $$\psi(q,r,t)$$
            % - Parameter centerXSpline: 1-D `TensorSpline` for $$x_c(t)$$
            % - Parameter centerYSpline: 1-D `TensorSpline` for $$y_c(t)$$
            % - Returns self: `TensorSplineStreamfunction` instance
            arguments
                streamfunctionSpline (1,1) TensorSpline
                centerXSpline (1,1) TensorSpline
                centerYSpline (1,1) TensorSpline
            end

            if streamfunctionSpline.numDimensions ~= 3
                error("TensorSplineStreamfunction:InvalidStreamfunctionSpline", ...
                    "streamfunctionSpline must be a 3-D TensorSpline with dimensions [q r t].");
            end
            if centerXSpline.numDimensions ~= 1 || centerYSpline.numDimensions ~= 1
                error("TensorSplineStreamfunction:InvalidCenterSpline", ...
                    "centerXSpline and centerYSpline must be 1-D TensorSpline objects.");
            end

            self.streamfunctionSpline = streamfunctionSpline;
            self.centerXSpline = centerXSpline;
            self.centerYSpline = centerYSpline;
            self.name = "Tensor-spline streamfunction";
        end

        function values = xCenter(self, t)
            % Evaluate the fitted center trajectory $$x_c(t)$$.
            %
            % - Topic: Inspect fitted components
            % - Declaration: values = xCenter(self,t)
            % - Parameter self: `TensorSplineStreamfunction` instance
            % - Parameter t: evaluation times in seconds
            % - Returns values: fitted center x-position in meters
            arguments
                self
                t {mustBeNumeric, mustBeReal}
            end

            values = self.centerXSpline.valueAtPoints(t);
        end

        function values = yCenter(self, t)
            % Evaluate the fitted center trajectory $$y_c(t)$$.
            %
            % - Topic: Inspect fitted components
            % - Declaration: values = yCenter(self,t)
            % - Parameter self: `TensorSplineStreamfunction` instance
            % - Parameter t: evaluation times in seconds
            % - Returns values: fitted center y-position in meters
            arguments
                self
                t {mustBeNumeric, mustBeReal}
            end

            values = self.centerYSpline.valueAtPoints(t);
        end

        function values = uBackground(self, t)
            % Evaluate the background x-velocity $$d x_c / dt$$.
            %
            % - Topic: Evaluate fitted diagnostics
            % - Declaration: values = uBackground(self,t)
            % - Parameter self: `TensorSplineStreamfunction` instance
            % - Parameter t: evaluation times in seconds
            % - Returns values: background x-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
            end

            values = self.centerXSpline.valueAtPoints(t, D=1);
        end

        function values = vBackground(self, t)
            % Evaluate the background y-velocity $$d y_c / dt$$.
            %
            % - Topic: Evaluate fitted diagnostics
            % - Declaration: values = vBackground(self,t)
            % - Parameter self: `TensorSplineStreamfunction` instance
            % - Parameter t: evaluation times in seconds
            % - Returns values: background y-velocity in $$m s^{-1}$$
            arguments
                self
                t {mustBeNumeric, mustBeReal}
            end

            values = self.centerYSpline.valueAtPoints(t, D=1);
        end

        function psiValue = psi(self, t, x, y)
            % Evaluate the centered-frame streamfunction on absolute coordinates.
            %
            % The implementation shifts the query points to
            % `q = x - xCenter(t)` and `r = y - yCenter(t)` before
            % evaluating the stored streamfunction spline.
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
            % Evaluate the mesoscale x-velocity $$-\psi_r$$.
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
            % Evaluate the mesoscale y-velocity $$\psi_q$$.
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
            % Evaluate the normal strain field $$\sigma_n = -2\psi_{qr}$$.
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
            % Evaluate the shear strain field $$\sigma_s = \psi_{qq} - \psi_{rr}$$.
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
            % Evaluate the relative-vorticity field $$\zeta = \psi_{qq} + \psi_{rr}$$.
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
        function self = fitFromTrajectories(xCell, yCell, tCell, options)
            % Fit the center trajectory and centered-frame streamfunction from drifter positions.
            %
            % Pass one numeric vector per drifter in each of `xCell`,
            % `yCell`, and `tCell`. The center trajectory splines are
            % estimated first from all x and y observations, then each
            % drifter trajectory is differentiated with a cubic-or-lower
            % interpolating spline to obtain centered velocities. The final
            % streamfunction fit uses only the derivative relations
            % `u = -psi_r` and `v = psi_q`.
            %
            % The centered-frame streamfunction is represented as
            %
            % $$
            % \psi(q,r,t) = \sum_{k=1}^{M_t} \sum_{j=1}^{M_r} \sum_{i=1}^{M_q}
            % c_{ijk} B_i(q) B_j(r) T_k(t),
            % $$
            %
            % where `B_i`, `B_j`, and `T_k` are the tensor-product B-spline
            % basis functions. Because the fit uses only `u = -psi_r` and
            % `v = psi_q`, the coefficients are unchanged if
            %
            % $$
            % \psi(q,r,t) \mapsto \psi(q,r,t) + a(t).
            % $$
            %
            % The implementation removes this gauge directly in
            % coefficient space. Let `c_k \in \mathbb{R}^{M_q M_r}` be the
            % spatial coefficient block for one time basis function and
            %
            % $$
            % P = \begin{bmatrix}
            % I_{M_q M_r - 1} \\
            % -\mathbf{1}^T
            % \end{bmatrix},
            % \qquad
            % c_k = P \alpha_k.
            % $$
            %
            % Then every retained spatial coefficient block satisfies
            % `\mathbf{1}^T c_k = 0`, which removes the spatially constant
            % mode associated with `a(t)`. Stacking all time blocks gives
            %
            % $$
            % c = Z \alpha, \qquad Z = I_{M_t} \otimes P,
            % $$
            %
            % and the fitted coefficients come from the unconstrained
            % reduced least-squares problem
            %
            % $$
            % \hat{\alpha} = \arg\min_{\alpha} \|R Z \alpha - U\|_2^2,
            % \qquad
            % \hat{c} = Z \hat{\alpha},
            % $$
            %
            % where `R` is the spline derivative design matrix and `U`
            % stacks the centered drifter velocities.
            %
            % When `psiKnotPoints` is omitted or empty, the fitter uses the
            % default degree `psiS = [2 2 0]` and builds padded terminated
            % knot vectors from the centered observations.
            %
            % - Topic: Fit the model
            % - Declaration: self = fitFromTrajectories(xCell,yCell,tCell,options)
            % - Parameter xCell: cell array of drifter x-position vectors in meters
            % - Parameter yCell: cell array of drifter y-position vectors in meters
            % - Parameter tCell: cell array of drifter time vectors in seconds
            % - Parameter options.psiKnotPoints: optional cell array `{qKnot, rKnot, tKnot}` for the streamfunction basis
            % - Parameter options.psiS: optional spline degree vector `[Sq Sr St]`, default `[2 2 0]`
            % - Returns self: fitted `TensorSplineStreamfunction`
            arguments
                xCell cell
                yCell cell
                tCell cell
                options.psiKnotPoints = []
                options.psiS (1,3) double {mustBeInteger, mustBeNonnegative} = [2 2 0]
            end

            xCell = reshape(xCell, [], 1);
            yCell = reshape(yCell, [], 1);
            tCell = reshape(tCell, [], 1);

            if ~(numel(xCell) == numel(yCell) && numel(yCell) == numel(tCell))
                error("TensorSplineStreamfunction:TrajectoryCellCountMismatch", ...
                    "xCell, yCell, and tCell must contain the same number of drifter records.");
            end

            psiS = reshape(options.psiS, 1, []);

            for iDrifter = 1:numel(tCell)
                validateattributes(xCell{iDrifter}, {'numeric'}, {'vector', 'real', 'finite'});
                validateattributes(yCell{iDrifter}, {'numeric'}, {'vector', 'real', 'finite'});
                validateattributes(tCell{iDrifter}, {'numeric'}, {'vector', 'real', 'finite'});

                xCell{iDrifter} = reshape(xCell{iDrifter}, [], 1);
                yCell{iDrifter} = reshape(yCell{iDrifter}, [], 1);
                tCell{iDrifter} = reshape(tCell{iDrifter}, [], 1);

                if ~(numel(xCell{iDrifter}) == numel(yCell{iDrifter}) && ...
                        numel(yCell{iDrifter}) == numel(tCell{iDrifter}))
                    error("TensorSplineStreamfunction:TrajectoryLengthMismatch", ...
                        "Each drifter record must have matching x, y, and t lengths.");
                end
                if numel(tCell{iDrifter}) < 2
                    error("TensorSplineStreamfunction:TrajectoryTooShort", ...
                        "Each drifter record must contain at least two samples.");
                end
                if any(diff(tCell{iDrifter}) <= 0)
                    error("TensorSplineStreamfunction:NonmonotonicTime", ...
                        "Each drifter time vector must be strictly increasing.");
                end
            end

            nDrifters = numel(tCell);
            allT = vertcat(tCell{:});
            allX = vertcat(xCell{:});
            allY = vertcat(yCell{:});

            pooledTimes = sort(allT);
            tCenterData = pooledTimes(1:nDrifters:end);
            if tCenterData(end) ~= pooledTimes(end)
                tCenterData = [tCenterData; pooledTimes(end)];
            end

            centerS = min(numel(tCenterData), 4) - 1;
            centerKnotPoints = BSpline.knotPointsForDataPoints(tCenterData, S=centerS, splineDOF=numel(tCenterData));
            centerBasisMatrix = BSpline.matrixForDataPoints(allT, knotPoints=centerKnotPoints, S=centerS);
            centerXSpline = TensorSpline(S=centerS, knotPoints=centerKnotPoints, xi=centerBasisMatrix \ allX);
            centerYSpline = TensorSpline(S=centerS, knotPoints=centerKnotPoints, xi=centerBasisMatrix \ allY);

            numObservations = sum(cellfun(@numel, tCell));
            qAll = zeros(numObservations, 1);
            rAll = zeros(numObservations, 1);
            tAll = zeros(numObservations, 1);
            qdotAll = zeros(numObservations, 1);
            rdotAll = zeros(numObservations, 1);

            offset = 0;
            for iDrifter = 1:nDrifters
                ti = tCell{iDrifter};
                xi = xCell{iDrifter};
                yi = yCell{iDrifter};
                drifterS = min(numel(ti), 4) - 1;

                xSpline = InterpolatingSpline(ti, xi, S=drifterS);
                ySpline = InterpolatingSpline(ti, yi, S=drifterS);

                xc = centerXSpline.valueAtPoints(ti);
                yc = centerYSpline.valueAtPoints(ti);
                ub = centerXSpline.valueAtPoints(ti, D=1);
                vb = centerYSpline.valueAtPoints(ti, D=1);
                dxdt = xSpline.valueAtPoints(ti, D=1);
                dydt = ySpline.valueAtPoints(ti, D=1);

                indices = offset + (1:numel(ti));
                qAll(indices) = xi - xc;
                rAll(indices) = yi - yc;
                tAll(indices) = ti;
                qdotAll(indices) = dxdt - ub;
                rdotAll(indices) = dydt - vb;
                offset = offset + numel(ti);
            end

            if isempty(options.psiKnotPoints)
                psiKnotPoints = TensorSplineStreamfunction.defaultPsiKnotPoints(qAll, rAll, tAll, psiS);
            else
                psiKnotPoints = options.psiKnotPoints;
                if ~iscell(psiKnotPoints) || numel(psiKnotPoints) ~= 3
                    error("TensorSplineStreamfunction:InvalidPsiKnotPoints", ...
                        "psiKnotPoints must be a cell array {qKnot, rKnot, tKnot}.");
                end

                psiKnotPoints = reshape(psiKnotPoints, 1, []);
                for iDim = 1:3
                    validateattributes(psiKnotPoints{iDim}, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
                    psiKnotPoints{iDim} = reshape(psiKnotPoints{iDim}, [], 1);
                    if any(diff(psiKnotPoints{iDim}) < 0)
                        error("TensorSplineStreamfunction:InvalidPsiKnotPoints", ...
                            "Each psi knot vector must be nondecreasing.");
                    end
                end
            end

            try
                psiBasisSize = TensorSpline.basisSizeFromKnotCell(psiKnotPoints, psiS + 1);
            catch
                error("TensorSplineStreamfunction:InvalidPsiKnotPoints", ...
                    "psiKnotPoints must define at least one basis function in each dimension.");
            end

            qDomain = [psiKnotPoints{1}(1), psiKnotPoints{1}(end)];
            rDomain = [psiKnotPoints{2}(1), psiKnotPoints{2}(end)];
            tDomain = [psiKnotPoints{3}(1), psiKnotPoints{3}(end)];

            if any(qAll < qDomain(1) | qAll > qDomain(2))
                error("TensorSplineStreamfunction:ObservationOutsidePsiDomain", ...
                    "Centered q observations must lie inside the supplied psi spline domain.");
            end
            if any(rAll < rDomain(1) | rAll > rDomain(2))
                error("TensorSplineStreamfunction:ObservationOutsidePsiDomain", ...
                    "Centered r observations must lie inside the supplied psi spline domain.");
            end
            if any(tAll < tDomain(1) | tAll > tDomain(2))
                error("TensorSplineStreamfunction:ObservationOutsidePsiDomain", ...
                    "Observation times must lie inside the supplied psi time domain.");
            end

            pointMatrix = [qAll, rAll, tAll];
            Ru = -sparse(TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[0 1 0]));
            Rv = sparse(TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[1 0 0]));
            R = [Ru; Rv];
            U = [qdotAll; rdotAll];
            numSpatialCoefficients = prod(psiBasisSize(1:2));
            numTimeCoefficients = psiBasisSize(3);
            if numSpatialCoefficients == 1
                reducedSpatialBasis = sparse(1, 0);
                reducedBasisMatrix = sparse(prod(psiBasisSize), 0);
                reducedDesignMatrix = sparse(size(R, 1), 0);
                xi = zeros(prod(psiBasisSize), 1);
            else
                reducedSpatialBasis = spalloc(numSpatialCoefficients, numSpatialCoefficients - 1, ...
                    2 * (numSpatialCoefficients - 1));
                reducedSpatialBasis(1:(numSpatialCoefficients - 1), :) = speye(numSpatialCoefficients - 1);
                reducedSpatialBasis(numSpatialCoefficients, :) = -1;
                reducedBasisMatrix = kron(speye(numTimeCoefficients), reducedSpatialBasis);
                reducedDesignMatrix = R * reducedBasisMatrix;
                reducedCoefficients = reducedDesignMatrix \ U;
                xi = reducedBasisMatrix * reducedCoefficients;
            end

            streamfunctionSpline = TensorSpline(S=psiS, knotPoints=psiKnotPoints, xi=xi);

            self = TensorSplineStreamfunction(streamfunctionSpline, centerXSpline, centerYSpline);
            self.fitState = struct( ...
                "centerSupportTimes", tCenterData, ...
                "centerDegree", centerS, ...
                "centerKnotPoints", centerKnotPoints, ...
                "centerBasisMatrix", centerBasisMatrix, ...
                "psiS", psiS, ...
                "psiKnotPoints", psiKnotPoints, ...
                "psiBasisSize", psiBasisSize, ...
                "psiReducedSpatialBasis", reducedSpatialBasis, ...
                "psiReducedBasisMatrix", reducedBasisMatrix, ...
                "psiReducedDesignMatrix", reducedDesignMatrix, ...
                "psiDesignMatrix", R, ...
                "psiObservedVelocities", U);
        end

        function self = fitFromSynchronousTrajectories(x, y, t, options)
            % Fit the model from synchronous gridded drifter trajectories.
            %
            % Use this convenience wrapper when all drifters are observed at
            % the same times, stored as `x=[nT nDrifters]`,
            % `y=[nT nDrifters]`, and `t=[nT 1]`. The wrapper converts each
            % drifter column into the canonical cell-array input expected by
            % `fitFromTrajectories(...)`.
            %
            % - Topic: Fit the model
            % - Declaration: self = fitFromSynchronousTrajectories(x,y,t,options)
            % - Parameter x: synchronous drifter x positions in meters
            % - Parameter y: synchronous drifter y positions in meters
            % - Parameter t: shared observation-time vector in seconds
            % - Parameter options.psiKnotPoints: optional cell array `{qKnot, rKnot, tKnot}` for the streamfunction basis
            % - Parameter options.psiS: optional spline degree vector `[Sq Sr St]`, default `[2 2 0]`
            % - Returns self: fitted `TensorSplineStreamfunction`
            arguments
                x (:,:) double {mustBeReal, mustBeFinite}
                y (:,:) double {mustBeReal, mustBeFinite}
                t (:,1) double {mustBeReal, mustBeFinite}
                options.psiKnotPoints = []
                options.psiS (1,3) double {mustBeInteger, mustBeNonnegative} = [2 2 0]
            end

            if ~isequal(size(x), size(y))
                error("TensorSplineStreamfunction:SynchronousSizeMismatch", ...
                    "x and y must have the same size for synchronous trajectories.");
            end
            if size(x, 1) ~= numel(t)
                error("TensorSplineStreamfunction:SynchronousTimeSizeMismatch", ...
                    "t must have one row per synchronous observation time.");
            end

            xCell = cell(size(x, 2), 1);
            yCell = cell(size(y, 2), 1);
            tCell = cell(size(x, 2), 1);
            for iDrifter = 1:size(x, 2)
                xCell{iDrifter} = x(:, iDrifter);
                yCell{iDrifter} = y(:, iDrifter);
                tCell{iDrifter} = t;
            end

            self = TensorSplineStreamfunction.fitFromTrajectories(xCell, yCell, tCell, ...
                psiKnotPoints=options.psiKnotPoints, psiS=options.psiS);
        end
    end

    methods (Access = private)
        function [tEval, q, r] = centeredCoordinates(self, t, x, y)
            if ~isequal(size(x), size(y))
                error("TensorSplineStreamfunction:EvaluationSizeMismatch", ...
                    "x and y must have the same size.");
            end

            if isscalar(t)
                tEval = repmat(t, size(x));
            elseif isvector(t) && ismatrix(x) && size(x, 1) == numel(t)
                tEval = repmat(reshape(t, [], 1), 1, size(x, 2));
            else
                if ~isequal(size(t), size(x))
                    error("TensorSplineStreamfunction:EvaluationTimeSizeMismatch", ...
                        "t must be scalar or have the same size as x and y.");
                end
                tEval = t;
            end

            q = x - self.xCenter(tEval);
            r = y - self.yCenter(tEval);
        end
    end

    methods (Static, Access = private)
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
