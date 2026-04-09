classdef GriddedStreamfunctionUnitTests < matlab.unittest.TestCase

    methods (Test)
        function constructorFitsAndReturnsHandle(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            fit = GriddedStreamfunction(trajectories);

            testCase.verifyClass(fit, "GriddedStreamfunction")
            testCase.verifyTrue(isa(fit, "handle"))
            testCase.verifyNotEmpty(fit.streamfunctionSpline)
            testCase.verifyNotEmpty(fit.observedTrajectories)
            testCase.verifyNotEmpty(fit.centerOfMassTrajectory)
            testCase.verifyNotEmpty(fit.backgroundTrajectory)
            testCase.verifyTrue(isstruct(fit.decomposition))
            testCase.verifyTrue(isfield(fit.decomposition, "fixedFrame"))
            testCase.verifyTrue(isfield(fit.decomposition, "centeredFrame"))
            testCase.verifyEqual(numel(fit.observedTrajectories), numel(trajectories))
            testCase.verifyEqual(numel(fit.decomposition.fixedFrame.background), numel(trajectories))
            testCase.verifyEqual(numel(fit.decomposition.fixedFrame.mesoscale), numel(trajectories))
            testCase.verifyEqual(numel(fit.decomposition.fixedFrame.submesoscale), numel(trajectories))
            testCase.verifyEqual(numel(fit.decomposition.centeredFrame.mesoscale), numel(trajectories))
            testCase.verifyEqual(numel(fit.decomposition.centeredFrame.submesoscale), numel(trajectories))
            testCase.verifyEqual(fit.fastS, 3)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.x(t) - mean(x, 2))), 2e-1)
            testCase.verifyLessThanOrEqual(max(abs(fit.centerOfMassTrajectory.y(t) - mean(y, 2))), 2e-1)

            for iTrajectory = 1:numel(trajectories)
                testCase.verifyEqual(fit.observedTrajectories(iTrajectory).t, trajectories(iTrajectory).t)
            end
        end

        function removedPropertiesAreNotStored(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();

            fit = GriddedStreamfunction(trajectories);
            removedProperties = [ ...
                "sampleCounts", "observationTimes", "observedX", "observedY", ...
                "observedXVelocity", "observedYVelocity", "centeredX", "centeredY", ...
                "centerVelocityX", "centerVelocityY", "uMesoscaleObserved", ...
                "vMesoscaleObserved", "uMesoscaleComObserved", "vMesoscaleComObserved", ...
                "uMesoscaleRelativeObserved", "vMesoscaleRelativeObserved", ...
                "uBackgroundObserved", "vBackgroundObserved", "uSubmesoscaleObserved", ...
                "vSubmesoscaleObserved", "backgroundTrajectories", "mesoscaleTrajectories", ...
                "submesoscaleTrajectories", "centeredMesoscaleTrajectories", ...
                "centeredSubmesoscaleTrajectories", "fitDiagnostics"];

            for iProperty = 1:numel(removedProperties)
                testCase.verifyFalse(isprop(fit, removedProperties(iProperty)))
            end
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

        function anchoredTrajectoryComponentsMatchConventionsAndReconstructPositions(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            for iTrajectory = 1:numel(trajectories)
                trajectory = fit.observedTrajectories(iTrajectory);
                ti = trajectory.t;
                observedX = trajectory.x(ti);
                observedY = trajectory.y(ti);
                [~, centeredX, centeredY] = fit.centeredCoordinates(ti, observedX, observedY);

                backgroundX = fit.decomposition.fixedFrame.background(iTrajectory).x(ti);
                backgroundY = fit.decomposition.fixedFrame.background(iTrajectory).y(ti);
                mesoscaleX = fit.decomposition.fixedFrame.mesoscale(iTrajectory).x(ti);
                mesoscaleY = fit.decomposition.fixedFrame.mesoscale(iTrajectory).y(ti);
                submesoscaleX = fit.decomposition.fixedFrame.submesoscale(iTrajectory).x(ti);
                submesoscaleY = fit.decomposition.fixedFrame.submesoscale(iTrajectory).y(ti);
                centeredMesoscaleX = fit.decomposition.centeredFrame.mesoscale(iTrajectory).x(ti);
                centeredMesoscaleY = fit.decomposition.centeredFrame.mesoscale(iTrajectory).y(ti);
                centeredSubmesoscaleX = fit.decomposition.centeredFrame.submesoscale(iTrajectory).x(ti);
                centeredSubmesoscaleY = fit.decomposition.centeredFrame.submesoscale(iTrajectory).y(ti);
                backgroundShiftedX = fit.backgroundTrajectory.x(ti) - fit.backgroundTrajectory.x(ti(1));
                backgroundShiftedY = fit.backgroundTrajectory.y(ti) - fit.backgroundTrajectory.y(ti(1));

                testCase.verifyEqual(backgroundX(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(backgroundY(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(submesoscaleX(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(submesoscaleY(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(centeredSubmesoscaleX(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(centeredSubmesoscaleY(1), 0, "AbsTol", 1e-12)
                testCase.verifyEqual(mesoscaleX(1), observedX(1), "AbsTol", 1e-12)
                testCase.verifyEqual(mesoscaleY(1), observedY(1), "AbsTol", 1e-12)
                testCase.verifyEqual(centeredMesoscaleX(1), centeredX(1), "AbsTol", 1e-12)
                testCase.verifyEqual(centeredMesoscaleY(1), centeredY(1), "AbsTol", 1e-12)
                testCase.verifyEqual(backgroundX, backgroundShiftedX, "AbsTol", 1e-12)
                testCase.verifyEqual(backgroundY, backgroundShiftedY, "AbsTol", 1e-12)

                testCase.verifyEqual(backgroundX + mesoscaleX + submesoscaleX, observedX, "AbsTol", 1e-2)
                testCase.verifyEqual(backgroundY + mesoscaleY + submesoscaleY, observedY, "AbsTol", 1e-2)
                testCase.verifyEqual(centeredMesoscaleX + centeredSubmesoscaleX, centeredX, "AbsTol", 1e-2)
                testCase.verifyEqual(centeredMesoscaleY + centeredSubmesoscaleY, centeredY, "AbsTol", 1e-2)
            end
        end

        function componentTrajectoryDerivativesMatchVelocityDecomposition(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            for iTrajectory = 1:numel(trajectories)
                trajectory = fit.observedTrajectories(iTrajectory);
                ti = trajectory.t;
                observedX = trajectory.x(ti);
                observedY = trajectory.y(ti);
                observedU = trajectory.u(ti);
                observedV = trajectory.v(ti);

                background = fit.decomposition.fixedFrame.background(iTrajectory);
                mesoscale = fit.decomposition.fixedFrame.mesoscale(iTrajectory);
                submesoscale = fit.decomposition.fixedFrame.submesoscale(iTrajectory);
                centeredMesoscale = fit.decomposition.centeredFrame.mesoscale(iTrajectory);
                centeredSubmesoscale = fit.decomposition.centeredFrame.submesoscale(iTrajectory);

                testCase.verifyEqual(background.u(ti), fit.uBackground(ti), "AbsTol", 1e-12)
                testCase.verifyEqual(background.v(ti), fit.vBackground(ti), "AbsTol", 1e-12)
                testCase.verifyEqual(mesoscale.u(ti), fit.uMesoscale(ti, observedX, observedY), "AbsTol", 1e-12)
                testCase.verifyEqual(mesoscale.v(ti), fit.vMesoscale(ti, observedX, observedY), "AbsTol", 1e-12)
                testCase.verifyEqual(background.u(ti) + mesoscale.u(ti) + submesoscale.u(ti), observedU, "AbsTol", 1e-12)
                testCase.verifyEqual(background.v(ti) + mesoscale.v(ti) + submesoscale.v(ti), observedV, "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    centeredMesoscale.u(ti) + centeredSubmesoscale.u(ti), ...
                    observedU - fit.centerOfMassTrajectory.u(ti), ...
                    "AbsTol", 1e-12)
                testCase.verifyEqual( ...
                    centeredMesoscale.v(ti) + centeredSubmesoscale.v(ti), ...
                    observedV - fit.centerOfMassTrajectory.v(ti), ...
                    "AbsTol", 1e-12)
            end
        end

        function constructorReadsTrajectoryEvaluationsInsteadOfDataValues(testCase)
            trajectories = GriddedStreamfunctionUnitTests.smoothedTrajectoryData();
            fit = GriddedStreamfunction(trajectories, psiKnotPoints=GriddedStreamfunctionUnitTests.zeroPsiKnotPoints());

            maxObservedRawDifference = 0;
            for iTrajectory = 1:numel(trajectories)
                trajectory = fit.observedTrajectories(iTrajectory);
                ti = trajectory.t;
                observedX = trajectory.x(ti);
                observedY = trajectory.y(ti);
                rawX = trajectory.x.dataValues;

                reconstructedX = fit.decomposition.fixedFrame.background(iTrajectory).x(ti) + ...
                    fit.decomposition.fixedFrame.mesoscale(iTrajectory).x(ti) + ...
                    fit.decomposition.fixedFrame.submesoscale(iTrajectory).x(ti);
                reconstructedY = fit.decomposition.fixedFrame.background(iTrajectory).y(ti) + ...
                    fit.decomposition.fixedFrame.mesoscale(iTrajectory).y(ti) + ...
                    fit.decomposition.fixedFrame.submesoscale(iTrajectory).y(ti);

                testCase.verifyEqual(reconstructedX, observedX, "AbsTol", 1e-2)
                testCase.verifyEqual(reconstructedY, observedY, "AbsTol", 1e-2)
                maxObservedRawDifference = max(maxObservedRawDifference, max(abs(observedX - rawX)));
            end

            testCase.verifyGreaterThan(maxObservedRawDifference, 1e-3)
        end

        function reducedBasisGaugeRemovesScalarSpatialMode(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            basisSize = fit.streamfunctionSpline.basisSize;
            xiByTime = reshape(fit.streamfunctionSpline.xi, prod(basisSize(1:2)), basisSize(3));
            gaugeModeMatrix = GriddedStreamfunctionUnitTests.constantGaugeMode(fit.psiKnotPoints, fit.psiS, basisSize);

            testCase.verifyEqual(size(gaugeModeMatrix, 2), 1)
            testCase.verifyEqual(gaugeModeMatrix.' * xiByTime, zeros(1, basisSize(3)), "AbsTol", 1e-10)
        end

        function evaluationInputShapes(testCase)
            [~, t, x, y, trajectories] = GriddedStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = GriddedStreamfunction(trajectories);

            tGrid = repmat(t, 1, size(x, 2));
            mxGrid = repmat(fit.centerOfMassTrajectory.x(t), 1, size(x, 2));
            myGrid = repmat(fit.centerOfMassTrajectory.y(t), 1, size(x, 2));
            scalarTime = t(4);
            scalarRow = repmat(scalarTime, 1, size(x, 2));

            [tEval, q, r] = fit.centeredCoordinates(t, x, y);

            testCase.verifyEqual(tEval, tGrid, "AbsTol", 1e-12)
            testCase.verifyEqual(q, x - mxGrid, "AbsTol", 1e-12)
            testCase.verifyEqual(r, y - myGrid, "AbsTol", 1e-12)
            testCase.verifyEqual(fit.uMesoscale(scalarTime, x(4, :), y(4, :)), fit.uMesoscale(scalarRow, x(4, :), y(4, :)), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.psiMesoscale(t, x, y), fit.psiMesoscale(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.uMesoscale(t, x, y), fit.uMesoscale(tGrid, x, y), "AbsTol", 1e-12)
            testCase.verifyEqual(fit.vMesoscale(t, x, y), fit.vMesoscale(tGrid, x, y), "AbsTol", 1e-12)

            testCase.verifyError(@() fit.centeredCoordinates(t, x(1:(end - 1), :), y), "GriddedStreamfunction:EvaluationSizeMismatch")
            testCase.verifyError(@() fit.centeredCoordinates(t(1:(end - 1)), x, y), "GriddedStreamfunction:EvaluationTimeSizeMismatch")
            testCase.verifyError(@() fit.uMesoscale(t, x(1:(end - 1), :), y), "GriddedStreamfunction:EvaluationSizeMismatch")
            testCase.verifyError(@() fit.uMesoscale(t(1:(end - 1)), x, y), "GriddedStreamfunction:EvaluationTimeSizeMismatch")
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

        function gaugeModeMatrix = constantGaugeMode(psiKnotPoints, psiS, basisSize)
            numSpatialCoefficients = prod(basisSize(1:2));
            if numSpatialCoefficients == 1
                gaugeModeMatrix = 1;
                return
            end

            qDomain = [psiKnotPoints{1}(1), psiKnotPoints{1}(end)];
            rDomain = [psiKnotPoints{2}(1), psiKnotPoints{2}(end)];
            spatialSupportCounts = max(2 * basisSize(1:2), basisSize(1:2) + 2);
            qGaugePoints = linspace(qDomain(1), qDomain(2), spatialSupportCounts(1)).';
            rGaugePoints = linspace(rDomain(1), rDomain(2), spatialSupportCounts(2)).';
            [qGaugeGrid, rGaugeGrid] = ndgrid(qGaugePoints, rGaugePoints);
            spatialPointMatrix = [qGaugeGrid(:), rGaugeGrid(:)];
            spatialBasisMatrix = TensorSpline.matrixForPointMatrix(spatialPointMatrix, knotPoints=psiKnotPoints(1:2), S=psiS(1:2));
            constantGaugeVector = spatialBasisMatrix \ ones(size(qGaugeGrid(:)));
            gaugeModeMatrix = reshape(constantGaugeVector, [], 1);
        end
    end
end
