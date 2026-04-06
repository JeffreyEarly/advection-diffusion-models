classdef GriddedStreamfunctionUnitTests < matlab.unittest.TestCase

    methods (Test)
        function constructorFitsAndReturnsHandle(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            fit = GriddedStreamfunction(trajectories);

            testCase.verifyClass(fit, "GriddedStreamfunction")
            testCase.verifyTrue(isa(fit, "handle"))
            testCase.verifyNotEmpty(fit.streamfunctionSpline)
            testCase.verifyNotEmpty(fit.centerOfMassTrajectory)
            testCase.verifyNotEmpty(fit.backgroundVelocityTrajectory)
            testCase.verifyEqual(numel(fit.backgroundTrajectories), numel(trajectories))
            testCase.verifyEqual(numel(fit.mesoscaleTrajectories), numel(trajectories))
            testCase.verifyEqual(numel(fit.submesoscaleTrajectories), numel(trajectories))
            testCase.verifyEqual(numel(fit.centeredMesoscaleTrajectories), numel(trajectories))
            testCase.verifyEqual(numel(fit.centeredSubmesoscaleTrajectories), numel(trajectories))
            testCase.verifyEqual(fit.fastS, 3)
            testCase.verifyEqual(numel(fit.observationTimes), numel(t) * size(x, 2))
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.x(t) - mean(x, 2))), 2e-1)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.y(t) - mean(y, 2))), 2e-1)
        end

        function constructorRejectsRawSampleInputs(testCase)
            x = randn(10, 3);

            testCase.verifyError(@() GriddedStreamfunction(x), "MATLAB:validators:mustBeA")
        end

        function constructorRejectsCellArrayInputs(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            cellInput = cell(numel(trajectories), 1);
            for iTrajectory = 1:numel(trajectories)
                cellInput{iTrajectory} = trajectories(iTrajectory);
            end

            testCase.verifyError(@() GriddedStreamfunction(cellInput), "MATLAB:validators:mustBeA")
        end

        function synchronousLinearFieldRecovery(testCase)
            [model, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            fit = GriddedStreamfunction(trajectories);
            tGrid = repmat(t, 1, size(x, 2));
            uBackgroundGrid = repmat(fit.uBackground(t), 1, size(x, 2));
            vBackgroundGrid = repmat(fit.vBackground(t), 1, size(x, 2));
            uResolved = fit.uMesoscale(tGrid, x, y) + uBackgroundGrid;
            vResolved = fit.vMesoscale(tGrid, x, y) + vBackgroundGrid;

            testCase.verifyEqual(fit.representativeTimes, t, "AbsTol", 1e-12)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.x(t) - mean(x, 2))), 2e-1)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.y(t) - mean(y, 2))), 2e-1)
            testCase.verifyLessThanOrEqual(max(abs(fit.sigma_n(tGrid, x, y) - model.sigma_n), [], "all"), 4e-7)
            testCase.verifyLessThanOrEqual(max(abs(fit.sigma_s(tGrid, x, y) - model.sigma_s), [], "all"), 3e-7)
            testCase.verifyLessThanOrEqual(max(abs(fit.zeta(tGrid, x, y) - model.zeta), [], "all"), 2e-7)
            testCase.verifyLessThanOrEqual(max(abs(uResolved - model.u(tGrid, x, y)), [], "all"), 3e-4)
            testCase.verifyLessThanOrEqual(max(abs(vResolved - model.v(tGrid, x, y)), [], "all"), 3e-4)
        end

        function asynchronousLinearFieldRecovery(testCase)
            [model, trajectories] = GriddedStreamfunctionUnitTests.asynchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            uResolvedError = zeros(numel(trajectories), 1);
            vResolvedError = zeros(numel(trajectories), 1);
            sigmaNError = zeros(numel(trajectories), 1);
            sigmaSError = zeros(numel(trajectories), 1);
            zetaError = zeros(numel(trajectories), 1);

            for iTrajectory = 1:numel(trajectories)
                trajectory = trajectories(iTrajectory);
                ti = trajectory.t;
                xi = trajectory.x(ti);
                yi = trajectory.y(ti);

                uResolvedError(iTrajectory) = max(abs(fit.uMesoscale(ti, xi, yi) + fit.uBackground(ti) - model.u(ti, xi, yi)));
                vResolvedError(iTrajectory) = max(abs(fit.vMesoscale(ti, xi, yi) + fit.vBackground(ti) - model.v(ti, xi, yi)));
                sigmaNError(iTrajectory) = max(abs(fit.sigma_n(ti, xi, yi) - model.sigma_n));
                sigmaSError(iTrajectory) = max(abs(fit.sigma_s(ti, xi, yi) - model.sigma_s));
                zetaError(iTrajectory) = max(abs(fit.zeta(ti, xi, yi) - model.zeta));
            end

            testCase.verifyLessThanOrEqual(max(uResolvedError), 5e-4)
            testCase.verifyLessThanOrEqual(max(vResolvedError), 5e-4)
            testCase.verifyLessThanOrEqual(max(sigmaNError), 2e-7)
            testCase.verifyLessThanOrEqual(max(sigmaSError), 2e-7)
            testCase.verifyLessThanOrEqual(max(zetaError), 2e-7)
        end

        function mesoscaleSolveMatchesStoredDesignAndCoefficients(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            tData = linspace(t(1), t(end), 5).';
            timeData = (tData(1:(end - 1)) + tData(2:end))/2;
            fastS = 3;
            fastKnotPoints = BSpline.knotPointsForDataPoints(timeData, S=fastS);
            fastKnotPoints(1:(fastS + 1)) = tData(1);
            fastKnotPoints((end - fastS):end) = tData(end);
            q = x - mean(x, 2);
            r = y - mean(y, 2);
            psiKnotPoints = { ...
                GriddedStreamfunctionUnitTests.paddedQuadraticDomain(q) ...
                GriddedStreamfunctionUnitTests.paddedQuadraticDomain(r) ...
                fastKnotPoints ...
                };
            fit = GriddedStreamfunction(trajectories, psiKnotPoints=psiKnotPoints, ...
                psiS=[2 2 fastS], fastKnotPoints=fastKnotPoints, fastS=fastS);

            residual = fit.fitDiagnostics.mesoscaleDesignMatrix * fit.fitDiagnostics.psiReducedCoefficients - ...
                fit.fitDiagnostics.mesoscaleRightHandSide;

            testCase.verifyLessThanOrEqual(norm(residual, inf), 5e-7)
            testCase.verifyEqual(size(fit.fitDiagnostics.mesoscaleDesignMatrix, 2), numel(fit.fitDiagnostics.psiReducedCoefficients))
        end

        function comVelocityMatchesStoredDerivativeSolve(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            testCase.verifyEqual(fit.centerVelocityX, fit.fitDiagnostics.fastDerivativeMatrix * fit.fitDiagnostics.centerXCoefficients, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.centerVelocityY, fit.fitDiagnostics.fastDerivativeMatrix * fit.fitDiagnostics.centerYCoefficients, "AbsTol", 1e-12)
        end

        function backgroundObservedSamplesMatchComResidual(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            testCase.verifyEqual(fit.uBackgroundObserved, fit.centerVelocityX - fit.uMesoscaleComObserved, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.vBackgroundObserved, fit.centerVelocityY - fit.vMesoscaleComObserved, "AbsTol", 1e-12)
        end

        function velocityDecompositionReconstructsObservedSamples(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            testCase.verifyEqual( ...
                fit.uMesoscaleObserved + fit.uBackgroundObserved + fit.uSubmesoscaleObserved, ...
                fit.observedXVelocity, "AbsTol", 1e-12)
            testCase.verifyEqual( ...
                fit.vMesoscaleObserved + fit.vBackgroundObserved + fit.vSubmesoscaleObserved, ...
                fit.observedYVelocity, "AbsTol", 1e-12)
            testCase.verifyEqual( ...
                fit.uMesoscaleComObserved + fit.uMesoscaleRelativeObserved, ...
                fit.uMesoscaleObserved, "AbsTol", 1e-12)
            testCase.verifyEqual( ...
                fit.vMesoscaleComObserved + fit.vMesoscaleRelativeObserved, ...
                fit.vMesoscaleObserved, "AbsTol", 1e-12)
        end

        function anchoredTrajectoryComponentsMatchConventionsAndReconstructPositions(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            sampleStartIndex = 1;
            for iTrajectory = 1:numel(trajectories)
                ti = trajectories(iTrajectory).t;
                nSamples = numel(ti);
                sampleIndices = sampleStartIndex:(sampleStartIndex + nSamples - 1);

                backgroundX = fit.backgroundTrajectories(iTrajectory).x(ti);
                backgroundY = fit.backgroundTrajectories(iTrajectory).y(ti);
                mesoscaleX = fit.mesoscaleTrajectories(iTrajectory).x(ti);
                mesoscaleY = fit.mesoscaleTrajectories(iTrajectory).y(ti);
                submesoscaleX = fit.submesoscaleTrajectories(iTrajectory).x(ti);
                submesoscaleY = fit.submesoscaleTrajectories(iTrajectory).y(ti);
                centeredMesoscaleX = fit.centeredMesoscaleTrajectories(iTrajectory).x(ti);
                centeredMesoscaleY = fit.centeredMesoscaleTrajectories(iTrajectory).y(ti);
                centeredSubmesoscaleX = fit.centeredSubmesoscaleTrajectories(iTrajectory).x(ti);
                centeredSubmesoscaleY = fit.centeredSubmesoscaleTrajectories(iTrajectory).y(ti);

                testCase.verifyEqual(backgroundX(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(backgroundY(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(submesoscaleX(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(submesoscaleY(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(centeredSubmesoscaleX(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(centeredSubmesoscaleY(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(mesoscaleX(1), fit.observedX(sampleIndices(1)), "AbsTol", 1e-12)
                testCase.verifyEqual(mesoscaleY(1), fit.observedY(sampleIndices(1)), "AbsTol", 1e-12)
                testCase.verifyEqual(centeredMesoscaleX(1), fit.centeredX(sampleIndices(1)), "AbsTol", 1e-12)
                testCase.verifyEqual(centeredMesoscaleY(1), fit.centeredY(sampleIndices(1)), "AbsTol", 1e-12)

                testCase.verifyEqual(backgroundX + mesoscaleX + submesoscaleX, fit.observedX(sampleIndices), "AbsTol", 1e-2)
                testCase.verifyEqual(backgroundY + mesoscaleY + submesoscaleY, fit.observedY(sampleIndices), "AbsTol", 1e-2)
                testCase.verifyEqual(centeredMesoscaleX + centeredSubmesoscaleX, fit.centeredX(sampleIndices), "AbsTol", 1e-2)
                testCase.verifyEqual(centeredMesoscaleY + centeredSubmesoscaleY, fit.centeredY(sampleIndices), "AbsTol", 1e-2)

                sampleStartIndex = sampleStartIndex + nSamples;
            end
        end

        function componentTrajectoryDerivativesMatchStoredVelocitySamples(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            sampleStartIndex = 1;
            for iTrajectory = 1:numel(trajectories)
                ti = trajectories(iTrajectory).t;
                nSamples = numel(ti);
                sampleIndices = sampleStartIndex:(sampleStartIndex + nSamples - 1);

                testCase.verifyEqual( ...
                    fit.backgroundTrajectories(iTrajectory).x.valueAtPoints(ti, D=1), ...
                    fit.uBackground(ti), "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    fit.backgroundTrajectories(iTrajectory).y.valueAtPoints(ti, D=1), ...
                    fit.vBackground(ti), "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    fit.mesoscaleTrajectories(iTrajectory).x.valueAtPoints(ti, D=1), ...
                    fit.uMesoscaleObserved(sampleIndices), "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    fit.mesoscaleTrajectories(iTrajectory).y.valueAtPoints(ti, D=1), ...
                    fit.vMesoscaleObserved(sampleIndices), "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    fit.submesoscaleTrajectories(iTrajectory).x.valueAtPoints(ti, D=1), ...
                    fit.uSubmesoscaleObserved(sampleIndices), "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    fit.submesoscaleTrajectories(iTrajectory).y.valueAtPoints(ti, D=1), ...
                    fit.vSubmesoscaleObserved(sampleIndices), "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    fit.centeredMesoscaleTrajectories(iTrajectory).x.valueAtPoints(ti, D=1), ...
                    fit.uMesoscaleRelativeObserved(sampleIndices), "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    fit.centeredMesoscaleTrajectories(iTrajectory).y.valueAtPoints(ti, D=1), ...
                    fit.vMesoscaleRelativeObserved(sampleIndices), "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    fit.centeredSubmesoscaleTrajectories(iTrajectory).x.valueAtPoints(ti, D=1), ...
                    fit.fitDiagnostics.uCenteredSubmesoscaleObserved(sampleIndices), "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    fit.centeredSubmesoscaleTrajectories(iTrajectory).y.valueAtPoints(ti, D=1), ...
                    fit.fitDiagnostics.vCenteredSubmesoscaleObserved(sampleIndices), "AbsTol", 1e-12)

                sampleStartIndex = sampleStartIndex + nSamples;
            end
        end

        function constructorReadsTrajectoryEvaluationsInsteadOfDataValues(testCase)
            trajectories = GriddedStreamfunctionUnitTests.smoothedTrajectoryData();
            fit = GriddedStreamfunction(trajectories, psiKnotPoints=GriddedStreamfunctionUnitTests.zeroPsiKnotPoints());

            expectedObservedX = zeros(0, 1);
            expectedObservedY = zeros(0, 1);
            rawX = zeros(0, 1);
            for iTrajectory = 1:numel(trajectories)
                expectedObservedX = [expectedObservedX; trajectories(iTrajectory).x(trajectories(iTrajectory).t)]; %#ok<AGROW>
                expectedObservedY = [expectedObservedY; trajectories(iTrajectory).y(trajectories(iTrajectory).t)]; %#ok<AGROW>
                rawX = [rawX; trajectories(iTrajectory).x.dataValues]; %#ok<AGROW>
            end

            testCase.verifyEqual(fit.observedX, expectedObservedX, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.observedY, expectedObservedY, "AbsTol", 1e-12)
            testCase.verifyGreaterThan(max(abs(expectedObservedX - rawX)), 1e-3)
        end

        function reducedBasisGaugeRemovesScalarSpatialMode(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            basisSize = fit.streamfunctionSpline.basisSize;
            xiByTime = reshape(fit.streamfunctionSpline.xi, prod(basisSize(1:2)), basisSize(3));
            gaugeModeMatrix = fit.fitDiagnostics.psiGaugeModeMatrix;

            testCase.verifyEqual(size(gaugeModeMatrix, 2), 1)
            testCase.verifyEqual(gaugeModeMatrix.' * xiByTime, zeros(1, basisSize(3)), "AbsTol", 1e-10)
        end

        function evaluationInputShapes(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            tGrid = repmat(t, 1, size(x, 2));
            scalarTime = t(4);
            scalarRow = repmat(scalarTime, 1, size(x, 2));

            testCase.verifyEqual(fit.uMesoscale(scalarTime, x(4, :), y(4, :)), fit.uMesoscale(scalarRow, x(4, :), y(4, :)), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.psiMesoscale(t, x, y), fit.psiMesoscale(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.uMesoscale(t, x, y), fit.uMesoscale(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.vMesoscale(t, x, y), fit.vMesoscale(tGrid, x, y), "AbsTol", 1e-12)

            testCase.verifyError(@() fit.uMesoscale(t, x(1:end-1, :), y), "GriddedStreamfunction:EvaluationSizeMismatch")
            testCase.verifyError(@() fit.uMesoscale(t(1:end-1), x, y), "GriddedStreamfunction:EvaluationTimeSizeMismatch")
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

        function trajectories = trajectorySplinesFromMatrices(t, x, y, splineDegree)
            if nargin < 4
                splineDegree = 3;
            end

            trajectories = TrajectorySpline.empty(0, 1);
            for iDrifter = 1:size(x, 2)
                trajectories(end + 1, 1) = TrajectorySpline(t, x(:, iDrifter), y(:, iDrifter), S=splineDegree); %#ok<AGROW>
            end
        end

        function knotPoints = paddedQuadraticDomain(values)
            extent = max(abs(values), [], "all");
            padding = max(50, 0.25 * max(extent, eps));
            lowerBound = -(extent + padding);
            upperBound = extent + padding;
            knotPoints = [repmat(lowerBound, 3, 1); repmat(upperBound, 3, 1)];
        end

        function trajectories = smoothedTrajectoryData()
            t = linspace(0, 1, 21).';
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
