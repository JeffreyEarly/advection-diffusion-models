classdef GriddedStreamfunctionUnitTests < matlab.unittest.TestCase

    methods (Test)
        function constructorStoresTrajectoryComponents(testCase)
            t = linspace(0, 1, 5)';
            knotPoints = {[-1; -1; 1; 1], [-1; -1; 1; 1], [0; 0; 1; 1]};
            streamfunctionSpline = TensorSpline(S=[1 1 1], knotPoints=knotPoints, xi=zeros(8, 1));
            centerTrajectory = TrajectorySpline(t, 2*t, -t, S=3);
            backgroundTrajectory = TrajectorySpline(t, 0.1 + 0*t, -0.2 + 0*t, S=3);

            fit = GriddedStreamfunction(streamfunctionSpline, centerTrajectory, backgroundTrajectory);

            testCase.verifyClass(fit, "GriddedStreamfunction")
            testCase.verifySameHandle(fit.streamfunctionSpline, streamfunctionSpline)
            testCase.verifySameHandle(fit.centerOfMassTrajectory, centerTrajectory)
            testCase.verifySameHandle(fit.backgroundVelocityTrajectory, backgroundTrajectory)
        end

        function fitRejectsRawSampleInputs(testCase)
            x = randn(10, 3);

            testCase.verifyError(@() GriddedStreamfunction.fitFromTrajectorySplines(x), ...
                "MATLAB:validators:mustBeA")
        end

        function fitRejectsCellArrayInputs(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            cellInput = cell(numel(trajectories), 1);
            for iTrajectory = 1:numel(trajectories)
                cellInput{iTrajectory} = trajectories(iTrajectory);
            end

            testCase.verifyError(@() GriddedStreamfunction.fitFromTrajectorySplines(cellInput), ...
                "MATLAB:validators:mustBeA")
        end

        function synchronousLinearFieldRecovery(testCase)
            [model, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            fit = GriddedStreamfunction.fitFromTrajectorySplines(trajectories);
            tGrid = repmat(t, 1, size(x, 2));

            testCase.verifyEqual(fit.fitState.representativeTimes, t, "AbsTol", 1e-12)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.x(t) - mean(x, 2))), 2e-1)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.y(t) - mean(y, 2))), 2e-1)
            testCase.verifyLessThanOrEqual(max(abs(fit.uBackground(t) - model.u0)), 2e-3)
            testCase.verifyLessThanOrEqual(max(abs(fit.vBackground(t) - model.v0)), 1e-2)
            testCase.verifyLessThanOrEqual(max(abs(fit.sigma_n(tGrid, x, y) - model.sigma_n), [], "all"), 2e-8)
            testCase.verifyLessThanOrEqual(max(abs(fit.sigma_s(tGrid, x, y) - model.sigma_s), [], "all"), 2e-8)
            testCase.verifyLessThanOrEqual(max(abs(fit.zeta(tGrid, x, y) - model.zeta), [], "all"), 2e-8)
            testCase.verifyLessThanOrEqual(max(abs(fit.u(tGrid, x, y) - model.u(tGrid, x, y)), [], "all"), 3e-4)
            testCase.verifyLessThanOrEqual(max(abs(fit.v(tGrid, x, y) - model.v(tGrid, x, y)), [], "all"), 3e-4)
        end

        function asynchronousLinearFieldRecovery(testCase)
            [model, trajectories] = GriddedStreamfunctionUnitTests.asynchronousLinearFieldData();
            fit = GriddedStreamfunction.fitFromTrajectorySplines(trajectories);

            uBackgroundError = zeros(numel(trajectories), 1);
            vBackgroundError = zeros(numel(trajectories), 1);
            uTotalError = zeros(numel(trajectories), 1);
            vTotalError = zeros(numel(trajectories), 1);
            sigmaNError = zeros(numel(trajectories), 1);
            sigmaSError = zeros(numel(trajectories), 1);
            zetaError = zeros(numel(trajectories), 1);

            for iTrajectory = 1:numel(trajectories)
                trajectory = trajectories(iTrajectory);
                ti = trajectory.t;
                xi = trajectory.x(ti);
                yi = trajectory.y(ti);

                uBackgroundError(iTrajectory) = max(abs(fit.uBackground(ti) - model.u0));
                vBackgroundError(iTrajectory) = max(abs(fit.vBackground(ti) - model.v0));
                uTotalError(iTrajectory) = max(abs(fit.u(ti, xi, yi) - model.u(ti, xi, yi)));
                vTotalError(iTrajectory) = max(abs(fit.v(ti, xi, yi) - model.v(ti, xi, yi)));
                sigmaNError(iTrajectory) = max(abs(fit.sigma_n(ti, xi, yi) - model.sigma_n));
                sigmaSError(iTrajectory) = max(abs(fit.sigma_s(ti, xi, yi) - model.sigma_s));
                zetaError(iTrajectory) = max(abs(fit.zeta(ti, xi, yi) - model.zeta));
            end

            testCase.verifyLessThanOrEqual(max(uBackgroundError), 2e-3)
            testCase.verifyLessThanOrEqual(max(vBackgroundError), 1e-2)
            testCase.verifyLessThanOrEqual(max(uTotalError), 5e-4)
            testCase.verifyLessThanOrEqual(max(vTotalError), 5e-4)
            testCase.verifyLessThanOrEqual(max(sigmaNError), 3e-8)
            testCase.verifyLessThanOrEqual(max(sigmaSError), 3e-8)
            testCase.verifyLessThanOrEqual(max(zetaError), 3e-8)
        end

        function backgroundMatchesFixedFrameResidualProjection(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fitFromTrajectorySplines(trajectories);
            tObserved = fit.fitState.observationTimes;

            expectedRhoX = fit.fitState.observedXVelocity - fit.fitState.uMesoscaleObserved;
            expectedRhoY = fit.fitState.observedYVelocity - fit.fitState.vMesoscaleObserved;
            expectedBackgroundX = fit.fitState.fastBasisMatrix * fit.fitState.backgroundXCoefficients;
            expectedBackgroundY = fit.fitState.fastBasisMatrix * fit.fitState.backgroundYCoefficients;

            testCase.verifyEqual(fit.fitState.rhoX, expectedRhoX, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.fitState.rhoY, expectedRhoY, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.fitState.uBackgroundObserved, expectedBackgroundX, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.fitState.vBackgroundObserved, expectedBackgroundY, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.uBackground(tObserved), fit.fitState.uBackgroundObserved, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.vBackground(tObserved), fit.fitState.vBackgroundObserved, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.fitState.submesoscaleX, fit.fitState.rhoX - fit.fitState.uBackgroundObserved, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.fitState.submesoscaleY, fit.fitState.rhoY - fit.fitState.vBackgroundObserved, "AbsTol", 1e-12)
        end

        function fitterReadsTrajectoryEvaluationsInsteadOfDataValues(testCase)
            trajectories = GriddedStreamfunctionUnitTests.smoothedTrajectoryData();
            fit = GriddedStreamfunction.fitFromTrajectorySplines(trajectories, psiKnotPoints=GriddedStreamfunctionUnitTests.zeroPsiKnotPoints());

            expectedObservedX = zeros(0, 1);
            expectedObservedY = zeros(0, 1);
            rawX = zeros(0, 1);
            for iTrajectory = 1:numel(trajectories)
                expectedObservedX = [expectedObservedX; trajectories(iTrajectory).x(trajectories(iTrajectory).t)]; %#ok<AGROW>
                expectedObservedY = [expectedObservedY; trajectories(iTrajectory).y(trajectories(iTrajectory).t)]; %#ok<AGROW>
                rawX = [rawX; trajectories(iTrajectory).x.dataValues]; %#ok<AGROW>
            end

            testCase.verifyEqual(fit.fitState.observedX, expectedObservedX, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.fitState.observedY, expectedObservedY, "AbsTol", 1e-12)
            testCase.verifyGreaterThan(max(abs(expectedObservedX - rawX)), 1e-3)
        end

        function reducedBasisGaugeHasZeroSpatialCoefficientSum(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fitFromTrajectorySplines(trajectories);

            basisSize = fit.streamfunctionSpline.basisSize;
            xiByTime = reshape(fit.streamfunctionSpline.xi, prod(basisSize(1:2)), basisSize(3));

            testCase.verifyEqual(sum(xiByTime, 1), zeros(1, basisSize(3)), "AbsTol", 1e-12)
        end

        function evaluationInputShapes(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction.fitFromTrajectorySplines(trajectories);

            tGrid = repmat(t, 1, size(x, 2));
            scalarTime = t(4);
            scalarRow = repmat(scalarTime, 1, size(x, 2));

            testCase.verifyEqual(fit.u(scalarTime, x(4, :), y(4, :)), fit.u(scalarRow, x(4, :), y(4, :)), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.psi(t, x, y), fit.psi(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.u(t, x, y), fit.u(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.v(t, x, y), fit.v(tGrid, x, y), "AbsTol", 1e-12)

            testCase.verifyError(@() fit.u(t, x(1:end-1, :), y), ...
                "GriddedStreamfunction:EvaluationSizeMismatch")
            testCase.verifyError(@() fit.u(t(1:end-1), x, y), ...
                "GriddedStreamfunction:EvaluationTimeSizeMismatch")
        end
    end

    methods (Static, Access = private)
        function [model, t, x, y, trajectories] = synchronousLinearFieldData()
            sigma = 4.0e-6;
            theta = pi/9;
            zeta = -1.5e-6;
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
            integrator = AdvectionDiffusionIntegrator(model, 0);

            [x0, y0] = ndgrid([-600 0 600] + 1200, [-400 0 400] - 800);
            [t, x, y] = integrator.particleTrajectories(x0(:), y0(:), 12*3600, 900);
            trajectories = GriddedStreamfunctionUnitTests.trajectorySplinesFromMatrices(t, x, y);
        end

        function [model, trajectories] = asynchronousLinearFieldData()
            sigma = 4.0e-6;
            theta = pi/9;
            zeta = -1.5e-6;
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
            integrator = AdvectionDiffusionIntegrator(model, 0);

            [x0, y0] = ndgrid([-600 0 600] + 1200, [-400 0 400] - 800);
            [tFine, xFine, yFine] = integrator.particleTrajectories(x0(:), y0(:), 12*3600, 300);

            nDrifters = size(xFine, 2);
            trajectories = TrajectorySpline.empty(nDrifters, 0);
            trajectories(nDrifters, 1) = TrajectorySpline();
            drifterOffset = [-1; 0; 1; 1; 0; 1; 1; 0; -1];
            baseIndices = 3:3:(numel(tFine) - 2);
            jitterPattern = [0; 1; 0; -1; 0; 1; 0; -1; 0; 0];
            jitter = repmat(jitterPattern, ceil(numel(baseIndices) / numel(jitterPattern)), 1);
            jitter = jitter(1:numel(baseIndices));

            for iDrifter = 1:nDrifters
                indices = baseIndices + transpose(jitter) + drifterOffset(iDrifter);
                trajectories(iDrifter) = TrajectorySpline(tFine(indices), xFine(indices, iDrifter), yFine(indices, iDrifter), S=3);
            end
        end

        function trajectories = trajectorySplinesFromMatrices(t, x, y, splineDegree)
            if nargin < 4
                splineDegree = 3;
            end

            nTrajectories = size(x, 2);
            trajectories = TrajectorySpline.empty(nTrajectories, 0);
            trajectories(nTrajectories, 1) = TrajectorySpline();
            for iDrifter = 1:size(x, 2)
                trajectories(iDrifter) = TrajectorySpline(t, x(:, iDrifter), y(:, iDrifter), S=splineDegree);
            end
        end

        function trajectories = smoothedTrajectoryData()
            t = linspace(0, 1, 21)';
            x1 = 2 + 0.3 * t;
            y1 = -1 + 0.2 * t;
            x2 = 3 + 0.3 * t;
            y2 = -0.5 + 0.2 * t;
            perturbation = 0.05 * sin(2*pi*t);

            xSpline1 = ConstrainedSpline(t, x1 + perturbation, S=3, splineDOF=6);
            ySpline1 = ConstrainedSpline(t, y1 - perturbation, S=3, splineDOF=6);
            xSpline2 = ConstrainedSpline(t, x2 - perturbation, S=3, splineDOF=6);
            ySpline2 = ConstrainedSpline(t, y2 + perturbation, S=3, splineDOF=6);

            trajectories = TrajectorySpline.empty(2, 0);
            trajectories(2, 1) = TrajectorySpline();
            trajectories(1) = TrajectorySpline.fromComponentSplines(t, xSpline1, ySpline1);
            trajectories(2) = TrajectorySpline.fromComponentSplines(t, xSpline2, ySpline2);
        end

        function psiKnotPoints = zeroPsiKnotPoints()
            psiKnotPoints = {[-1; -1; -1; 1; 1; 1], [-1; -1; -1; 1; 1; 1], [0; 1]};
        end
    end
end
