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
            fit = GriddedStreamfunction.fitFromTrajectorySplines(trajectories, ...
                psiKnotPoints=GriddedStreamfunctionUnitTests.zeroPsiKnotPoints());

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

        function legacyAlignedSynchronousLinearFieldMatchesEffectiveDiagnostics(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunctionUnitTests.legacyAlignedFit(t, x, y, trajectories, 1);
            diagnostics = GriddedStreamfunctionUnitTests.effectiveDiagnostics(fit, t, x, y);
            parameters = EstimateLinearVelocityFieldParameters(x, y, t, ...
                [ModelParameter.u0v0, ModelParameter.strain, ModelParameter.vorticity], 1);

            mx = fit.centerOfMassTrajectory.x(t);
            my = fit.centerOfMassTrajectory.y(t);

            testCase.verifyEqual(diagnostics.effective.uBackground, parameters.u0 + zeros(size(t)), "AbsTol", 2e-3)
            testCase.verifyEqual(diagnostics.effective.vBackground, parameters.v0 + zeros(size(t)), "AbsTol", 1e-2)
            testCase.verifyEqual(fit.sigma_n(t, mx, my), parameters.sigma_n + zeros(size(t)), "AbsTol", 1e-10)
            testCase.verifyEqual(fit.sigma_s(t, mx, my), parameters.sigma_s + zeros(size(t)), "AbsTol", 1e-10)
            testCase.verifyEqual(fit.zeta(t, mx, my), parameters.zeta + zeros(size(t)), "AbsTol", 1e-10)
        end

        function legacyAlignedTimeVaryingBackgroundMatchesEffectiveDiagnostics(testCase)
            [t, x, y, trajectories] = GriddedStreamfunctionUnitTests.timeVaryingBackgroundOnlyData();
            fit = GriddedStreamfunctionUnitTests.legacyAlignedFit(t, x, y, trajectories, 4);
            diagnostics = GriddedStreamfunctionUnitTests.effectiveDiagnostics(fit, t, x, y);
            parameters = EstimateLinearVelocityFieldParameters(x, y, t, ModelParameter.u0v0, 4);

            mx = fit.centerOfMassTrajectory.x(t);
            my = fit.centerOfMassTrajectory.y(t);

            testCase.verifyEqual(diagnostics.effective.uBackground, parameters.u0, "AbsTol", 1e-7)
            testCase.verifyEqual(diagnostics.effective.vBackground, parameters.v0, "AbsTol", 1e-7)
            testCase.verifyEqual(fit.sigma_n(t, mx, my), zeros(size(t)), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.sigma_s(t, mx, my), zeros(size(t)), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.zeta(t, mx, my), zeros(size(t)), "AbsTol", 1e-12)
        end

        function effectiveDiagnosticsZeroMeanMesoscaleAndReconstructVelocity(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunctionUnitTests.legacyAlignedFit(t, x, y, trajectories, 4);
            diagnostics = GriddedStreamfunctionUnitTests.effectiveDiagnostics(fit, t, x, y);

            testCase.verifyEqual(mean(diagnostics.effective.uMesoscale, 2), zeros(size(t)), "AbsTol", 1e-12)
            testCase.verifyEqual(mean(diagnostics.effective.vMesoscale, 2), zeros(size(t)), "AbsTol", 1e-12)
            testCase.verifyEqual( ...
                diagnostics.effective.uMesoscale + diagnostics.effective.uBackgroundObserved + diagnostics.effective.uSubmesoscale, ...
                diagnostics.raw.uObserved, "AbsTol", 1e-12)
            testCase.verifyEqual( ...
                diagnostics.effective.vMesoscale + diagnostics.effective.vBackgroundObserved + diagnostics.effective.vSubmesoscale, ...
                diagnostics.raw.vObserved, "AbsTol", 1e-12)
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
            trajectories = TrajectorySpline.empty(0, 1);
            drifterOffset = [-1; 0; 1; 1; 0; 1; 1; 0; -1];
            baseIndices = 3:3:(numel(tFine) - 2);
            jitterPattern = [0; 1; 0; -1; 0; 1; 0; -1; 0; 0];
            jitter = repmat(jitterPattern, ceil(numel(baseIndices) / numel(jitterPattern)), 1);
            jitter = jitter(1:numel(baseIndices));

            for iDrifter = 1:nDrifters
                indices = baseIndices + transpose(jitter) + drifterOffset(iDrifter);
                trajectories(end + 1, 1) = TrajectorySpline(tFine(indices), xFine(indices, iDrifter), yFine(indices, iDrifter), S=3); %#ok<AGROW>
            end
        end

        function [t, x, y, trajectories] = timeVaryingBackgroundOnlyData()
            t = linspace(0, 12*3600, 49).';
            T = t(end);

            xCenter = 0.05*t + 0.03*(t.^2)/(2*T);
            yCenter = -0.04*t + 0.02*(t.^2)/(2*T);

            [xOffset, yOffset] = ndgrid([-300 0 300], [-200 0 200]);
            x = xCenter + reshape(xOffset, 1, []);
            y = yCenter + reshape(yOffset, 1, []);
            trajectories = GriddedStreamfunctionUnitTests.trajectorySplinesFromMatrices(t, x, y);
        end

        function fit = legacyAlignedFit(t, x, y, trajectories, dof)
            options = GriddedStreamfunctionUnitTests.legacyAlignedFitOptions(t, x, y, dof);
            fit = GriddedStreamfunction.fitFromTrajectorySplines(trajectories, ...
                psiKnotPoints=options.psiKnotPoints, psiS=options.psiS, ...
                fastKnotPoints=options.fastKnotPoints, fastS=options.fastS);
        end

        function options = legacyAlignedFitOptions(t, x, y, dof)
            [timeBasisMatrix, tKnotPoints, timeSplineDegree] = GriddedStreamfunctionUnitTests.legacyTimeBasis(t, dof);
            psiS = [2 2 timeSplineDegree];

            q = x - timeBasisMatrix * (timeBasisMatrix \ mean(x, 2));
            r = y - timeBasisMatrix * (timeBasisMatrix \ mean(y, 2));

            options = struct( ...
                "timeBasisMatrix", timeBasisMatrix, ...
                "fastKnotPoints", tKnotPoints, ...
                "fastS", timeSplineDegree, ...
                "psiKnotPoints", {{ ...
                    GriddedStreamfunctionUnitTests.paddedQuadraticDomain(q) ...
                    GriddedStreamfunctionUnitTests.paddedQuadraticDomain(r) ...
                    tKnotPoints ...
                    }}, ...
                "psiS", psiS);
        end

        function [basisMatrix, knotPoints, splineDegree] = legacyTimeBasis(t, dof)
            if dof == 1
                splineDegree = 0;
                knotPoints = [t(1); t(end)];
                basisMatrix = ones(numel(t), 1);
                return
            end

            splineDegree = min([dof 4]) - 1;
            tData = linspace(t(1), t(end), dof + 1).';
            timeData = (tData(1:(end - 1)) + tData(2:end))/2;
            knotPoints = BSpline.knotPointsForDataPoints(timeData, S=splineDegree);
            knotPoints(1:(splineDegree + 1)) = tData(1);
            knotPoints((end - splineDegree):end) = tData(end);
            basisMatrix = BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=splineDegree);
        end

        function knotPoints = paddedQuadraticDomain(values)
            extent = max(abs(values), [], "all");
            padding = max(50, 0.25 * max(extent, eps));
            lowerBound = -(extent + padding);
            upperBound = extent + padding;
            knotPoints = [repmat(lowerBound, 3, 1); repmat(upperBound, 3, 1)];
        end

        function diagnostics = effectiveDiagnostics(fit, t, x, y)
            nTimes = numel(t);
            nDrifters = size(x, 2);

            uObserved = reshape(fit.fitState.observedXVelocity, nTimes, nDrifters);
            vObserved = reshape(fit.fitState.observedYVelocity, nTimes, nDrifters);
            uMesoscaleRaw = reshape(fit.fitState.uMesoscaleObserved, nTimes, nDrifters);
            vMesoscaleRaw = reshape(fit.fitState.vMesoscaleObserved, nTimes, nDrifters);
            uBackgroundRawObserved = reshape(fit.fitState.uBackgroundObserved, nTimes, nDrifters);
            vBackgroundRawObserved = reshape(fit.fitState.vBackgroundObserved, nTimes, nDrifters);
            uSubmesoscale = reshape(fit.fitState.submesoscaleX, nTimes, nDrifters);
            vSubmesoscale = reshape(fit.fitState.submesoscaleY, nTimes, nDrifters);

            uMesoscaleMean = mean(uMesoscaleRaw, 2);
            vMesoscaleMean = mean(vMesoscaleRaw, 2);

            diagnostics = struct( ...
                "raw", struct( ...
                    "uObserved", uObserved, ...
                    "vObserved", vObserved), ...
                "effective", struct( ...
                    "uMesoscale", uMesoscaleRaw - uMesoscaleMean, ...
                    "vMesoscale", vMesoscaleRaw - vMesoscaleMean, ...
                    "uBackground", fit.uBackground(t) + uMesoscaleMean, ...
                    "vBackground", fit.vBackground(t) + vMesoscaleMean, ...
                    "uBackgroundObserved", uBackgroundRawObserved + uMesoscaleMean, ...
                    "vBackgroundObserved", vBackgroundRawObserved + vMesoscaleMean, ...
                    "uSubmesoscale", uSubmesoscale, ...
                    "vSubmesoscale", vSubmesoscale));
        end

        function trajectories = trajectorySplinesFromMatrices(t, x, y, splineDegree)
            if nargin < 4
                splineDegree = 3;
            end

            trajectories = TrajectorySpline.empty(0, 1);
            for iDrifter = 1:size(x, 2)
                trajectories(end + 1, 1) = TrajectorySpline(t, x(:, iDrifter), y(:, iDrifter), S=splineDegree); %#ok<AGROW>
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

            trajectories = TrajectorySpline.empty(0, 1);
            trajectories(end + 1, 1) = TrajectorySpline.fromComponentSplines(t, xSpline1, ySpline1);
            trajectories(end + 1, 1) = TrajectorySpline.fromComponentSplines(t, xSpline2, ySpline2);
        end

        function psiKnotPoints = zeroPsiKnotPoints()
            psiKnotPoints = {[-1; -1; -1; 1; 1; 1], [-1; -1; -1; 1; 1; 1], [0; 1]};
        end
    end
end
