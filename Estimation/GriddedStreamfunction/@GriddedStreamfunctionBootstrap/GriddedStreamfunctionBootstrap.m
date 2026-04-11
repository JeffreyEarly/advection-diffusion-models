classdef GriddedStreamfunctionBootstrap < handle
    % Bootstrap whole-drifter gridded-streamfunction fits and consensus scores.
    %
    % `GriddedStreamfunctionBootstrap` fits one deterministic
    % `GriddedStreamfunction` on the full drifter ensemble, then constructs
    % a whole-drifter bootstrap ensemble by resampling the input
    % `TrajectorySpline` objects with replacement. Each replicate is
    % summarized by COM-local mesoscale diagnostics
    % $$u_c(t), v_c(t), \sigma_n(t), \sigma_s(t), \zeta(t)$$ evaluated at
    % the fitted center-of-mass trajectory, together with a scalar
    % diffusivity diagnostic
    % $$\kappa = \langle (x_{\mathrm{sm}}^2 + y_{\mathrm{sm}}^2)/(4\Delta t)\rangle.$$
    %
    % The bootstrap ensemble is ranked by a consensus score formed from
    % time-local `KernelDensityEstimate` fits to the bootstrap cloud. The
    % score uses a 2-D density fit for `(uCenter, vCenter)`, a 2-D
    % density fit for `(sigma_n, sigma_s)`, and a 1-D density fit for
    % `zeta`, omitting constrained blocks when `mesoscaleConstraint`
    % forces them to vanish.
    %
    % ```matlab
    % bootstrap = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=100);
    % bestFit = bootstrap.bestFit();
    % quantiles = bootstrap.summaryQuantiles([0.16 0.5 0.84]);
    % ```
    %
    % - Topic: Create a bootstrap ensemble
    % - Topic: Inspect bootstrap properties
    % - Topic: Reconstruct bootstrap fits
    % - Topic: Summarize bootstrap uncertainty
    % - Declaration: classdef GriddedStreamfunctionBootstrap < handle

    properties (SetAccess = private)
        % Full-data gridded-streamfunction fit used as the reference solution.
        %
        % `fullFit` is the deterministic `GriddedStreamfunction` fit on
        % the original input drifters before any resampling is applied.
        %
        % - Topic: Inspect bootstrap properties
        fullFit

        % Original drifter trajectories used to seed the bootstrap ensemble.
        %
        % `observedTrajectories` preserves the original
        % `TrajectorySpline` column vector in the user-supplied drifter
        % order.
        %
        % - Topic: Inspect bootstrap properties
        observedTrajectories

        % Whole-drifter resampling indices for each bootstrap replicate.
        %
        % `bootstrapIndices(i,:)` contains the drifter indices drawn with
        % replacement for bootstrap replicate `i`.
        %
        % - Topic: Inspect bootstrap properties
        bootstrapIndices

        % Times used to store COM-local bootstrap summaries.
        %
        % `queryTimes` is a column vector of times in seconds. When not
        % supplied explicitly, it is the sorted unique observation times
        % restricted to the common overlap interval of the original
        % drifters.
        %
        % - Topic: Inspect bootstrap properties
        queryTimes

        % Times used to compute the consensus score.
        %
        % `scoreTimes` is always an exact subset of `queryTimes`.
        %
        % - Topic: Inspect bootstrap properties
        scoreTimes

        % Full-data COM-local mesoscale summary evaluated on `queryTimes`.
        %
        % `fullSummary` contains the fields `uCenter`, `vCenter`,
        % `sigma_n`, `sigma_s`, and `zeta`, each stored as a column vector
        % aligned with `queryTimes`.
        %
        % - Topic: Inspect bootstrap properties
        fullSummary struct = struct()

        % Bootstrap COM-local mesoscale summaries evaluated on `queryTimes`.
        %
        % `summary.uCenter`, `summary.vCenter`, `summary.sigma_n`,
        % `summary.sigma_s`, and `summary.zeta` are arrays of size
        % `[numel(queryTimes) nBootstraps]`.
        %
        % - Topic: Inspect bootstrap properties
        summary struct = struct()

        % Full-data scalar bootstrap diagnostic summary.
        %
        % `fullScalarSummary.kappaEstimate` stores the scalar diffusivity
        % diagnostic computed from the full-data fitted submesoscale
        % trajectories.
        %
        % - Topic: Inspect bootstrap properties
        fullScalarSummary struct = struct()

        % Scalar diagnostic summaries for each bootstrap replicate.
        %
        % `scalarSummary.kappaEstimate` is a row vector of length
        % `nBootstraps`.
        %
        % - Topic: Inspect bootstrap properties
        scalarSummary struct = struct()

        % Consensus-score components and joint score for each bootstrap fit.
        %
        % `scores.uv`, `scores.strain`, `scores.zeta`, and `scores.joint`
        % are row vectors of length `nBootstraps`.
        %
        % - Topic: Inspect bootstrap properties
        scores struct = struct()

        % Exact reconstruction metadata for each bootstrap replicate.
        %
        % `bootstrapMetadata.fastKnotPoints{i}` stores the resolved fast
        % temporal knot vector used by bootstrap replicate `i`, and
        % `bootstrapMetadata.psiKnotPoints{i}` stores the corresponding
        % resolved mesoscale knot cell `{qKnot, rKnot, tKnot}`.
        %
        % - Topic: Reconstruct bootstrap fits
        bootstrapMetadata struct = struct()

        % Number of bootstrap replicates in the ensemble.
        %
        % - Topic: Inspect bootstrap properties
        nBootstraps (1,1) double {mustBeInteger, mustBeNonnegative} = 0

        % Random-number seed used for whole-drifter resampling.
        %
        % - Topic: Inspect bootstrap properties
        randomSeed (1,1) double {mustBeInteger, mustBeFinite} = 0
    end

    methods
        function self = GriddedStreamfunctionBootstrap(trajectories, options)
            % Create a whole-drifter bootstrap ensemble for `GriddedStreamfunction`.
            %
            % Use this constructor when the uncertainty analysis should
            % reflect sensitivity to which drifters were observed, rather
            % than only the residual variability within one fitted
            % ensemble.
            %
            % The constructor first fits the full-data
            % `GriddedStreamfunction`, then resamples the drifter
            % trajectories with replacement and fits one replicate per
            % bootstrap draw. Each replicate stores COM-local mesoscale
            % diagnostics on `queryTimes`, a scalar submesoscale
            % diffusivity diagnostic, and the resolved spline knot vectors
            % required to reconstruct the replicate exactly later.
            %
            % - Topic: Create a bootstrap ensemble
            % - Declaration: self = GriddedStreamfunctionBootstrap(trajectories,nBootstraps=...,randomSeed=...,queryTimes=...,scoreTimes=...,scoreStride=...,psiKnotPoints=...,psiS=...,fastKnotPoints=...,fastS=...,mesoscaleConstraint=...)
            % - Parameter trajectories: nonempty vector of `TrajectorySpline` drifters
            % - Parameter nBootstraps: number of whole-drifter bootstrap replicates, default `100`
            % - Parameter randomSeed: integer random seed used for resampling, default `0`
            % - Parameter queryTimes: optional strictly increasing query times for stored summaries
            % - Parameter scoreTimes: optional strictly increasing subset of `queryTimes` used for consensus scoring
            % - Parameter scoreStride: optional positive stride used to subsample `queryTimes` into `scoreTimes`, default `1`
            % - Parameter psiKnotPoints: optional cell array `{qKnot, rKnot, tKnot}` for the mesoscale basis
            % - Parameter psiS: optional mesoscale spline degree vector `[Sq Sr St]`, default `[2 2 0]`
            % - Parameter fastKnotPoints: optional fast temporal knot vector for COM and background
            % - Parameter fastS: optional fast temporal spline degree, default `3`
            % - Parameter mesoscaleConstraint: optional hard mesoscale constraint `"none"`, `"zeroVorticity"`, or `"zeroStrain"`
            % - Returns self: fitted bootstrap ensemble
            arguments (Input)
                trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
                options.nBootstraps (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.nBootstraps, 3)} = 100
                options.randomSeed (1,1) double {mustBeInteger, mustBeFinite} = 0
                options.queryTimes = []
                options.scoreTimes = []
                options.scoreStride (1,1) double {mustBeInteger, mustBePositive} = 1
                options.psiKnotPoints = []
                options.psiS (1,3) double {mustBeInteger, mustBeNonnegative} = [2 2 0]
                options.fastKnotPoints = []
                options.fastS (1,1) double {mustBeInteger, mustBeNonnegative} = 3
                options.mesoscaleConstraint {mustBeTextScalar, mustBeMember(options.mesoscaleConstraint, ["none", "zeroVorticity", "zeroStrain"])} = "none"
            end

            trajectories = reshape(trajectories, [], 1);
            fitOptions = struct( ...
                "psiKnotPoints", {options.psiKnotPoints}, ...
                "psiS", options.psiS, ...
                "fastKnotPoints", {options.fastKnotPoints}, ...
                "fastS", options.fastS, ...
                "mesoscaleConstraint", string(options.mesoscaleConstraint));
            fitArguments = namedargs2cell(fitOptions);
            fitOptionsWithoutDecomposition = fitOptions;
            fitOptionsWithoutDecomposition.buildDecomposition = false;
            fitArgumentsWithoutDecomposition = namedargs2cell(fitOptionsWithoutDecomposition);

            [queryTimes, scoreTimes, scoreIndices] = GriddedStreamfunctionBootstrap.resolveBootstrapTimes( ...
                trajectories, options.queryTimes, options.scoreTimes, options.scoreStride);

            self.observedTrajectories = trajectories;
            self.nBootstraps = options.nBootstraps;
            self.randomSeed = options.randomSeed;
            self.queryTimes = queryTimes;
            self.scoreTimes = scoreTimes;

            self.fullFit = GriddedStreamfunction(trajectories, fitArguments{:});
            [self.fullSummary, self.fullScalarSummary] = GriddedStreamfunctionBootstrap.extractSummaryFromFit( ...
                self.fullFit, queryTimes);

            nQuery = numel(queryTimes);
            nTrajectories = numel(trajectories);
            self.bootstrapIndices = zeros(self.nBootstraps, nTrajectories);
            self.summary = emptyTimeSummary(nQuery, self.nBootstraps);
            self.scalarSummary = emptyScalarSummary(self.nBootstraps);
            self.bootstrapMetadata = struct( ...
                "fastKnotPoints", {cell(self.nBootstraps, 1)}, ...
                "psiKnotPoints", {cell(self.nBootstraps, 1)});

            originalRng = rng;
            cleanupRng = onCleanup(@() rng(originalRng));
            rng(self.randomSeed);

            for iBootstrap = 1:self.nBootstraps
                sampledIndices = randi(nTrajectories, 1, nTrajectories);
                sampledTrajectories = reshape(trajectories(sampledIndices), [], 1);
                fit = GriddedStreamfunction(sampledTrajectories, fitArgumentsWithoutDecomposition{:});
                [summary, scalarSummary] = GriddedStreamfunctionBootstrap.extractSummaryFromFit(fit, queryTimes);

                self.bootstrapIndices(iBootstrap, :) = sampledIndices;
                self.summary.uCenter(:, iBootstrap) = summary.uCenter;
                self.summary.vCenter(:, iBootstrap) = summary.vCenter;
                self.summary.sigma_n(:, iBootstrap) = summary.sigma_n;
                self.summary.sigma_s(:, iBootstrap) = summary.sigma_s;
                self.summary.zeta(:, iBootstrap) = summary.zeta;
                self.scalarSummary.kappaEstimate(1, iBootstrap) = scalarSummary.kappaEstimate;
                self.bootstrapMetadata.fastKnotPoints{iBootstrap} = fit.fastKnotPoints;
                self.bootstrapMetadata.psiKnotPoints{iBootstrap} = fit.psiKnotPoints;
            end

            self.scores = GriddedStreamfunctionBootstrap.computeConsensusScores( ...
                self.summary, scoreIndices, self.fullFit.mesoscaleConstraint);
        end

        function index = bestBootstrapIndex(self)
            % Return the bootstrap index with the highest joint consensus score.
            %
            % - Topic: Reconstruct bootstrap fits
            % - Declaration: index = bestBootstrapIndex(self)
            % - Returns index: 1-based bootstrap index of the top-ranked replicate
            arguments (Input)
                self (1,1) GriddedStreamfunctionBootstrap
            end
            arguments (Output)
                index (1,1) double
            end

            [~, index] = max(self.scores.joint, [], 2);
        end

        function fit = fitForBootstrap(self, iBootstrap)
            % Reconstruct one bootstrap replicate exactly from saved metadata.
            %
            % Use this method to recover a full `GriddedStreamfunction`
            % object for a particular bootstrap replicate without storing
            % every replicate fit in memory.
            %
            % - Topic: Reconstruct bootstrap fits
            % - Declaration: fit = fitForBootstrap(self,iBootstrap)
            % - Parameter iBootstrap: 1-based bootstrap replicate index
            % - Returns fit: reconstructed `GriddedStreamfunction` replicate
            arguments (Input)
                self (1,1) GriddedStreamfunctionBootstrap
                iBootstrap (1,1) double {mustBeInteger, mustBePositive}
            end
            arguments (Output)
                fit (1,1) GriddedStreamfunction
            end

            if iBootstrap > self.nBootstraps
                error("GriddedStreamfunctionBootstrap:InvalidBootstrapIndex", ...
                    "iBootstrap must be between 1 and %d.", self.nBootstraps);
            end

            sampledTrajectories = reshape(self.observedTrajectories(self.bootstrapIndices(iBootstrap, :)), [], 1);
            psiKnotPoints = self.bootstrapMetadata.psiKnotPoints(iBootstrap);
            fastKnotPoints = self.bootstrapMetadata.fastKnotPoints(iBootstrap);
            fitOptions = struct( ...
                "psiKnotPoints", psiKnotPoints, ...
                "psiS", self.fullFit.psiS, ...
                "fastKnotPoints", fastKnotPoints, ...
                "fastS", self.fullFit.fastS, ...
                "mesoscaleConstraint", self.fullFit.mesoscaleConstraint);
            fitArguments = namedargs2cell(fitOptions);
            fit = GriddedStreamfunction(sampledTrajectories, fitArguments{:});
        end

        function fit = bestFit(self)
            % Reconstruct the top-ranked bootstrap replicate.
            %
            % - Topic: Reconstruct bootstrap fits
            % - Declaration: fit = bestFit(self)
            % - Returns fit: reconstructed `GriddedStreamfunction` for the best bootstrap replicate
            arguments (Input)
                self (1,1) GriddedStreamfunctionBootstrap
            end
            arguments (Output)
                fit (1,1) GriddedStreamfunction
            end

            fit = fitForBootstrap(self, bestBootstrapIndex(self));
        end

        function quantiles = summaryQuantiles(self, probabilities)
            % Return bootstrap quantiles for the stored summaries.
            %
            % The returned struct contains the time-varying summary fields
            % `uCenter`, `vCenter`, `sigma_n`, `sigma_s`, and `zeta` with
            % size `[numel(queryTimes) numel(probabilities)]`, together
            % with the scalar field `kappaEstimate` of size
            % `[1 numel(probabilities)]`.
            %
            % - Topic: Summarize bootstrap uncertainty
            % - Declaration: quantiles = summaryQuantiles(self,probabilities)
            % - Parameter probabilities: vector of quantile probabilities on `[0,1]`
            % - Returns quantiles: struct of bootstrap quantiles for the stored summary fields
            arguments (Input)
                self (1,1) GriddedStreamfunctionBootstrap
                probabilities {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
            end
            arguments (Output)
                quantiles struct
            end

            probabilities = reshape(probabilities, 1, []);
            if any(probabilities < 0 | probabilities > 1)
                error("GriddedStreamfunctionBootstrap:InvalidProbabilities", ...
                    "probabilities must lie on the closed interval [0, 1].");
            end

            quantiles = struct( ...
                "uCenter", quantile(self.summary.uCenter, probabilities, 2), ...
                "vCenter", quantile(self.summary.vCenter, probabilities, 2), ...
                "sigma_n", quantile(self.summary.sigma_n, probabilities, 2), ...
                "sigma_s", quantile(self.summary.sigma_s, probabilities, 2), ...
                "zeta", quantile(self.summary.zeta, probabilities, 2), ...
                "kappaEstimate", quantile(self.scalarSummary.kappaEstimate, probabilities, 2));
        end
    end

    methods (Static, Access = private)
        [queryTimes, scoreTimes, scoreIndices] = resolveBootstrapTimes(trajectories, queryTimesOption, scoreTimesOption, scoreStride)
        [summary, scalarSummary] = extractSummaryFromFit(fit, queryTimes)
        kappaEstimate = kappaEstimateFromFit(fit)
        scores = computeConsensusScores(summary, scoreIndices, mesoscaleConstraint)
    end
end

function summary = emptyTimeSummary(nQuery, nBootstraps)
summary = struct( ...
    "uCenter", zeros(nQuery, nBootstraps), ...
    "vCenter", zeros(nQuery, nBootstraps), ...
    "sigma_n", zeros(nQuery, nBootstraps), ...
    "sigma_s", zeros(nQuery, nBootstraps), ...
    "zeta", zeros(nQuery, nBootstraps));
end

function scalarSummary = emptyScalarSummary(nBootstraps)
scalarSummary = struct("kappaEstimate", zeros(1, nBootstraps));
end
