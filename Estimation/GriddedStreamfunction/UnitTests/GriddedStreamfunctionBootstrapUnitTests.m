classdef GriddedStreamfunctionBootstrapUnitTests < matlab.unittest.TestCase

    methods (Test)
        function constructorStoresOutputsAndReproducesSeed(testCase)
            [~, t, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrapA = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=6, randomSeed=17);
            bootstrapB = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=6, randomSeed=17);

            testCase.verifyClass(bootstrapA, "GriddedStreamfunctionBootstrap")
            testCase.verifyTrue(isa(bootstrapA, "handle"))
            testCase.verifyClass(bootstrapA.fullFit, "GriddedStreamfunction")
            testCase.verifyEqual(bootstrapA.nBootstraps, 6)
            testCase.verifyEqual(bootstrapA.randomSeed, 17)
            testCase.verifyEqual(size(bootstrapA.bootstrapIndices), [6 numel(trajectories)])
            testCase.verifyEqual(bootstrapA.queryTimes, t, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.scoreTimes, t, "AbsTol", 1e-12)
            testCase.verifySize(bootstrapA.summary.uCenter, [numel(t), 6])
            testCase.verifySize(bootstrapA.summary.sigma_n, [numel(t), 6])
            testCase.verifySize(bootstrapA.scalarSummary.kappaEstimate, [1 6])
            testCase.verifyTrue(all(isfinite(bootstrapA.scores.uv), "all"))
            testCase.verifyTrue(all(isfinite(bootstrapA.scores.strain), "all"))
            testCase.verifyTrue(all(isfinite(bootstrapA.scores.zeta), "all"))
            testCase.verifyTrue(all(isfinite(bootstrapA.scores.joint), "all"))
            testCase.verifyGreaterThanOrEqual(min(bootstrapA.bootstrapIndices, [], "all"), 1)
            testCase.verifyLessThanOrEqual(max(bootstrapA.bootstrapIndices, [], "all"), numel(trajectories))

            testCase.verifyEqual(bootstrapA.bootstrapIndices, bootstrapB.bootstrapIndices)
            testCase.verifyEqual(bootstrapA.summary.uCenter, bootstrapB.summary.uCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.summary.vCenter, bootstrapB.summary.vCenter, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.summary.sigma_n, bootstrapB.summary.sigma_n, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.summary.sigma_s, bootstrapB.summary.sigma_s, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.summary.zeta, bootstrapB.summary.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.scalarSummary.kappaEstimate, bootstrapB.scalarSummary.kappaEstimate, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.scores.uv, bootstrapB.scores.uv, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.scores.strain, bootstrapB.scores.strain, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.scores.zeta, bootstrapB.scores.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrapA.scores.joint, bootstrapB.scores.joint, "AbsTol", 1e-12)

            for iBootstrap = 1:bootstrapA.nBootstraps
                testCase.verifyEqual(bootstrapA.bootstrapMetadata.fastKnotPoints{iBootstrap}, ...
                    bootstrapB.bootstrapMetadata.fastKnotPoints{iBootstrap}, "AbsTol", 1e-12)
                for iDim = 1:3
                    testCase.verifyEqual(bootstrapA.bootstrapMetadata.psiKnotPoints{iBootstrap}{iDim}, ...
                        bootstrapB.bootstrapMetadata.psiKnotPoints{iBootstrap}{iDim}, "AbsTol", 1e-12)
                end
            end
        end

        function fullSummaryMatchesSynchronousLinearFieldTruth(testCase)
            [model, t, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrap = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=4, randomSeed=5);

            testCase.verifyLessThanOrEqual(max(abs(bootstrap.fullSummary.uCenter - model.u0)), 2e-3)
            testCase.verifyLessThanOrEqual(max(abs(bootstrap.fullSummary.vCenter - model.v0)), 5e-3)
            testCase.verifyLessThanOrEqual(max(abs(bootstrap.fullSummary.sigma_n - model.sigma_n)), 1e-12)
            testCase.verifyLessThanOrEqual(max(abs(bootstrap.fullSummary.sigma_s - model.sigma_s)), 1e-12)
            testCase.verifyLessThanOrEqual(max(abs(bootstrap.fullSummary.zeta - model.zeta)), 1e-12)
            testCase.verifyEqual(bootstrap.queryTimes, t, "AbsTol", 1e-12)
        end

        function bestFitAndStoredReplicateReconstructExactly(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrap = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=6, randomSeed=9);
            iBest = bootstrap.bestBootstrapIndex();
            fit = bootstrap.bestFit();
            [summary, scalarSummary] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit(fit, bootstrap.queryTimes);

            testCase.verifyEqual(fit.fastKnotPoints, bootstrap.bootstrapMetadata.fastKnotPoints{iBest}, "AbsTol", 1e-12)
            for iDim = 1:3
                testCase.verifyEqual(fit.psiKnotPoints{iDim}, bootstrap.bootstrapMetadata.psiKnotPoints{iBest}{iDim}, "AbsTol", 1e-12)
            end
            testCase.verifyEqual(summary.uCenter, bootstrap.summary.uCenter(:, iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.vCenter, bootstrap.summary.vCenter(:, iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.sigma_n, bootstrap.summary.sigma_n(:, iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.sigma_s, bootstrap.summary.sigma_s(:, iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.zeta, bootstrap.summary.zeta(:, iBest), "AbsTol", 1e-12)
            testCase.verifyEqual(scalarSummary.kappaEstimate, bootstrap.scalarSummary.kappaEstimate(iBest), "AbsTol", 1e-12)
        end

        function storedOutputsMatchFullyReconstructedFitsForAllReplicates(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrap = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=4, randomSeed=11);

            for iBootstrap = 1:bootstrap.nBootstraps
                fit = bootstrap.fitForBootstrap(iBootstrap);
                [summary, scalarSummary] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit(fit, bootstrap.queryTimes);

                testCase.verifyEqual(fit.fastKnotPoints, bootstrap.bootstrapMetadata.fastKnotPoints{iBootstrap}, "AbsTol", 1e-12)
                for iDim = 1:3
                    testCase.verifyEqual(fit.psiKnotPoints{iDim}, ...
                        bootstrap.bootstrapMetadata.psiKnotPoints{iBootstrap}{iDim}, "AbsTol", 1e-12)
                end

                testCase.verifyEqual(summary.uCenter, bootstrap.summary.uCenter(:, iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(summary.vCenter, bootstrap.summary.vCenter(:, iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(summary.sigma_n, bootstrap.summary.sigma_n(:, iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(summary.sigma_s, bootstrap.summary.sigma_s(:, iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(summary.zeta, bootstrap.summary.zeta(:, iBootstrap), "AbsTol", 1e-12)
                testCase.verifyEqual(scalarSummary.kappaEstimate, bootstrap.scalarSummary.kappaEstimate(iBootstrap), "AbsTol", 1e-12)
            end
        end

        function explicitKnotsAreStoredAndReusedDuringReconstruction(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            referenceFit = GriddedStreamfunction(trajectories);

            bootstrap = GriddedStreamfunctionBootstrap( ...
                trajectories, ...
                nBootstraps=4, ...
                randomSeed=4, ...
                fastKnotPoints=referenceFit.fastKnotPoints, ...
                psiKnotPoints=referenceFit.psiKnotPoints, ...
                psiS=referenceFit.psiS, ...
                fastS=referenceFit.fastS);

            for iBootstrap = 1:bootstrap.nBootstraps
                testCase.verifyEqual(bootstrap.bootstrapMetadata.fastKnotPoints{iBootstrap}, referenceFit.fastKnotPoints, "AbsTol", 1e-12)
                for iDim = 1:3
                    testCase.verifyEqual(bootstrap.bootstrapMetadata.psiKnotPoints{iBootstrap}{iDim}, ...
                        referenceFit.psiKnotPoints{iDim}, "AbsTol", 1e-12)
                end
            end

            fit = bootstrap.fitForBootstrap(2);
            testCase.verifyEqual(fit.fastKnotPoints, referenceFit.fastKnotPoints, "AbsTol", 1e-12)
            for iDim = 1:3
                testCase.verifyEqual(fit.psiKnotPoints{iDim}, referenceFit.psiKnotPoints{iDim}, "AbsTol", 1e-12)
            end
        end

        function defaultKnotsStillReconstructExactly(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrap = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=5, randomSeed=12);
            fit = bootstrap.fitForBootstrap(2);
            [summary, scalarSummary] = GriddedStreamfunctionBootstrapUnitTests.summaryFromFit(fit, bootstrap.queryTimes);

            testCase.verifyEqual(fit.fastKnotPoints, bootstrap.bootstrapMetadata.fastKnotPoints{2}, "AbsTol", 1e-12)
            for iDim = 1:3
                testCase.verifyEqual(fit.psiKnotPoints{iDim}, bootstrap.bootstrapMetadata.psiKnotPoints{2}{iDim}, "AbsTol", 1e-12)
            end
            testCase.verifyEqual(summary.uCenter, bootstrap.summary.uCenter(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.vCenter, bootstrap.summary.vCenter(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.sigma_n, bootstrap.summary.sigma_n(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.sigma_s, bootstrap.summary.sigma_s(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(summary.zeta, bootstrap.summary.zeta(:, 2), "AbsTol", 1e-12)
            testCase.verifyEqual(scalarSummary.kappaEstimate, bootstrap.scalarSummary.kappaEstimate(2), "AbsTol", 1e-12)
        end

        function defaultQueryTimesUseOriginalCommonOverlap(testCase)
            [trajectories, expectedQueryTimes, overlapBounds] = ...
                GriddedStreamfunctionBootstrapUnitTests.partialSupportTrajectoryData();

            bootstrap = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=4, randomSeed=3);

            testCase.verifyEqual(bootstrap.queryTimes, expectedQueryTimes, "AbsTol", 1e-12)
            testCase.verifyGreaterThanOrEqual(min(bootstrap.queryTimes), overlapBounds(1))
            testCase.verifyLessThanOrEqual(max(bootstrap.queryTimes), overlapBounds(2))
            testCase.verifyTrue(all(isfinite(bootstrap.summary.uCenter), "all"))
            testCase.verifyTrue(all(isfinite(bootstrap.summary.sigma_n), "all"))
        end

        function constrainedBootstrapSkipsDegenerateScoreBlocks(testCase)
            [trajectoriesVorticity, trajectoriesStrain] = ...
                GriddedStreamfunctionBootstrapUnitTests.constrainedTrajectoryData();

            zeroVorticity = GriddedStreamfunctionBootstrap( ...
                trajectoriesVorticity, nBootstraps=4, randomSeed=8, mesoscaleConstraint="zeroVorticity");
            zeroStrain = GriddedStreamfunctionBootstrap( ...
                trajectoriesStrain, nBootstraps=4, randomSeed=8, mesoscaleConstraint="zeroStrain");

            testCase.verifyEqual(zeroVorticity.scores.zeta, zeros(1, zeroVorticity.nBootstraps), "AbsTol", 1e-12)
            testCase.verifyEqual(zeroVorticity.scores.joint, ...
                zeroVorticity.scores.uv + zeroVorticity.scores.strain, "AbsTol", 1e-12)
            testCase.verifyTrue(all(isfinite(zeroVorticity.scores.joint)))

            testCase.verifyEqual(zeroStrain.scores.strain, zeros(1, zeroStrain.nBootstraps), "AbsTol", 1e-12)
            testCase.verifyEqual(zeroStrain.scores.joint, ...
                zeroStrain.scores.uv + zeroStrain.scores.zeta, "AbsTol", 1e-12)
            testCase.verifyTrue(all(isfinite(zeroStrain.scores.joint)))
        end

        function consensusScoresMatchDirectGaussianReference(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=4, randomSeed=5);
            scoreIndices = arrayfun(@(ti) find(bootstrap.queryTimes == ti, 1), bootstrap.scoreTimes);

            expected = GriddedStreamfunctionBootstrapUnitTests.directConsensusScores( ...
                bootstrap.summary, scoreIndices, bootstrap.fullFit.mesoscaleConstraint);

            testCase.verifyEqual(bootstrap.scores.uv, expected.uv, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.scores.strain, expected.strain, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.scores.zeta, expected.zeta, "AbsTol", 1e-12)
            testCase.verifyEqual(bootstrap.scores.joint, expected.joint, "AbsTol", 1e-12)
        end

        function summaryQuantilesReturnExpectedShapesAndMonotoneBands(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

            bootstrap = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=5, randomSeed=14);
            quantiles = bootstrap.summaryQuantiles([0.1 0.5 0.9]);

            testCase.verifySize(quantiles.uCenter, [numel(bootstrap.queryTimes), 3])
            testCase.verifySize(quantiles.vCenter, [numel(bootstrap.queryTimes), 3])
            testCase.verifySize(quantiles.sigma_n, [numel(bootstrap.queryTimes), 3])
            testCase.verifySize(quantiles.sigma_s, [numel(bootstrap.queryTimes), 3])
            testCase.verifySize(quantiles.zeta, [numel(bootstrap.queryTimes), 3])
            testCase.verifySize(quantiles.kappaEstimate, [1, 3])

            testCase.verifyTrue(all(quantiles.uCenter(:, 1) <= quantiles.uCenter(:, 2)))
            testCase.verifyTrue(all(quantiles.uCenter(:, 2) <= quantiles.uCenter(:, 3)))
            testCase.verifyTrue(all(quantiles.sigma_n(:, 1) <= quantiles.sigma_n(:, 2)))
            testCase.verifyTrue(all(quantiles.sigma_n(:, 2) <= quantiles.sigma_n(:, 3)))
            testCase.verifyTrue(all(quantiles.kappaEstimate(1) <= quantiles.kappaEstimate(2)))
            testCase.verifyTrue(all(quantiles.kappaEstimate(2) <= quantiles.kappaEstimate(3)))
        end

        function scalarDiffusivityDiagnosticIsFiniteAndIncreases(testCase)
            trajectoriesLow = GriddedStreamfunctionBootstrapUnitTests.diffusiveTrajectoryData(0.25, 21);
            trajectoriesHigh = GriddedStreamfunctionBootstrapUnitTests.diffusiveTrajectoryData(2.0, 21);

            bootstrapLow = GriddedStreamfunctionBootstrap(trajectoriesLow, nBootstraps=4, randomSeed=5);
            bootstrapHigh = GriddedStreamfunctionBootstrap(trajectoriesHigh, nBootstraps=4, randomSeed=5);

            testCase.verifyTrue(isfinite(bootstrapLow.fullScalarSummary.kappaEstimate))
            testCase.verifyTrue(isfinite(bootstrapHigh.fullScalarSummary.kappaEstimate))
            testCase.verifyGreaterThanOrEqual(bootstrapLow.fullScalarSummary.kappaEstimate, 0)
            testCase.verifyGreaterThanOrEqual(bootstrapHigh.fullScalarSummary.kappaEstimate, 0)
            testCase.verifyGreaterThan(bootstrapHigh.fullScalarSummary.kappaEstimate, ...
                bootstrapLow.fullScalarSummary.kappaEstimate)
            testCase.verifyGreaterThan(mean(bootstrapHigh.scalarSummary.kappaEstimate), ...
                mean(bootstrapLow.scalarSummary.kappaEstimate))
        end

        function scalarDiffusivityDiagnosticMatchesLegacySplineIntegration(testCase)
            [~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();
            bootstrap = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=4, randomSeed=5);

            legacyKappa = GriddedStreamfunctionBootstrapUnitTests.legacyKappaFromFit(bootstrap.fullFit);

            testCase.verifyEqual(bootstrap.fullScalarSummary.kappaEstimate, legacyKappa, "AbsTol", 1e-12)
        end
    end

    methods (Static)
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
            trajectories = GriddedStreamfunctionBootstrapUnitTests.trajectorySplinesFromMatrices(t, x, y);
        end

        function [trajectories, expectedQueryTimes, overlapBounds] = partialSupportTrajectoryData()
            sigma = 3.0e-6;
            theta = pi/10;
            zeta = -1.0e-6;
            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=0.04, v0=-0.02);
            integrator = AdvectionDiffusionIntegrator(model, 0);
            [x0, y0] = ndgrid([-400 0 400] + 800, [-300 0 300] - 500);
            [tFine, xFine, yFine] = integrator.particleTrajectories(x0(:), y0(:), 10*3600, 300);

            startIndices = [1; 4; 2; 6; 3];
            endIndices = [numel(tFine) - 2; numel(tFine); numel(tFine) - 4; numel(tFine) - 1; numel(tFine) - 3];
            nDrifters = numel(startIndices);
            trajectories = TrajectorySpline.empty(0, 1);
            allTimes = cell(nDrifters, 1);

            for iDrifter = 1:nDrifters
                indices = startIndices(iDrifter):3:endIndices(iDrifter);
                allTimes{iDrifter} = tFine(indices);
                trajectories(end + 1, 1) = TrajectorySpline.fromData(tFine(indices), xFine(indices, iDrifter), yFine(indices, iDrifter), S=3); %#ok<AGROW>
            end

            overlapBounds = [max(cellfun(@(ti) ti(1), allTimes)), min(cellfun(@(ti) ti(end), allTimes))];
            pooledTimes = unique(vertcat(allTimes{:}), "sorted");
            expectedQueryTimes = pooledTimes(pooledTimes >= overlapBounds(1) & pooledTimes <= overlapBounds(2));
        end

        function [trajectoriesZeroVorticity, trajectoriesZeroStrain] = constrainedTrajectoryData()
            integratorVorticity = AdvectionDiffusionIntegrator(LinearVelocityField(sigma=4.0e-6, theta=pi/8, zeta=0, u0=0.06, v0=-0.03), 0);
            integratorStrain = AdvectionDiffusionIntegrator(LinearVelocityField(sigma=0, theta=pi/8, zeta=-2.0e-6, u0=0.06, v0=-0.03), 0);
            [x0, y0] = ndgrid([-600 0 600] + 1200, [-400 0 400] - 800);
            [t, xVorticity, yVorticity] = integratorVorticity.particleTrajectories(x0(:), y0(:), 12*3600, 900);
            [~, xStrain, yStrain] = integratorStrain.particleTrajectories(x0(:), y0(:), 12*3600, 900);

            trajectoriesZeroVorticity = GriddedStreamfunctionBootstrapUnitTests.trajectorySplinesFromMatrices(t, xVorticity, yVorticity);
            trajectoriesZeroStrain = GriddedStreamfunctionBootstrapUnitTests.trajectorySplinesFromMatrices(t, xStrain, yStrain);
        end

        function trajectories = diffusiveTrajectoryData(kappa, rngSeed)
            sigma = 3.0e-6;
            theta = pi/12;
            zeta = -1.0e-6;
            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=0.03, v0=-0.02);
            integrator = AdvectionDiffusionIntegrator(model, kappa);
            integrator.stepSize = 900;

            [x0, y0] = ndgrid([-500 0 500] + 1000, [-350 0 350] - 600);
            rng(rngSeed);
            [t, x, y] = integrator.particleTrajectories(x0(:), y0(:), 18*3600, 1800);
            trajectories = GriddedStreamfunctionBootstrapUnitTests.trajectorySplinesFromMatrices(t, x, y);
        end

        function trajectories = trajectorySplinesFromMatrices(t, x, y, splineDegree)
            if nargin < 4
                splineDegree = 3;
            end

            trajectories = TrajectorySpline.empty(0, 1);
            for iDrifter = 1:size(x, 2)
                trajectories(end + 1, 1) = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=splineDegree); %#ok<AGROW>
            end
        end

        function [summary, scalarSummary] = summaryFromFit(fit, queryTimes)
            mx = reshape(fit.centerOfMassTrajectory.x(queryTimes), [], 1);
            my = reshape(fit.centerOfMassTrajectory.y(queryTimes), [], 1);

            summary = struct( ...
                "uCenter", reshape(fit.uMesoscale(queryTimes, mx, my), [], 1), ...
                "vCenter", reshape(fit.vMesoscale(queryTimes, mx, my), [], 1), ...
                "sigma_n", reshape(fit.sigma_n(queryTimes, mx, my), [], 1), ...
                "sigma_s", reshape(fit.sigma_s(queryTimes, mx, my), [], 1), ...
                "zeta", reshape(fit.zeta(queryTimes, mx, my), [], 1));
            scalarSummary = struct("kappaEstimate", GriddedStreamfunctionBootstrapUnitTests.kappaFromFit(fit));
        end

        function kappaEstimate = kappaFromFit(fit)
            nTrajectories = numel(fit.decomposition.fixedFrame.submesoscale);
            kappaValues = zeros(nTrajectories, 1);

            for iTrajectory = 1:nTrajectories
                trajectory = fit.decomposition.fixedFrame.submesoscale(iTrajectory);
                ti = trajectory.t;
                duration = ti(end) - ti(1);
                kappaValues(iTrajectory) = (trajectory.x(ti(end)).^2 + trajectory.y(ti(end)).^2) / (4 * duration);
            end

            kappaEstimate = mean(kappaValues);
        end

        function kappaEstimate = legacyKappaFromFit(fit)
            trajectories = reshape(fit.observedTrajectories, [], 1);
            nTrajectories = numel(trajectories);
            tCell = cell(nTrajectories, 1);
            xCell = cell(nTrajectories, 1);
            yCell = cell(nTrajectories, 1);
            xDotCell = cell(nTrajectories, 1);
            yDotCell = cell(nTrajectories, 1);

            for iTrajectory = 1:nTrajectories
                ti = reshape(trajectories(iTrajectory).t, [], 1);
                tCell{iTrajectory} = ti;
                xCell{iTrajectory} = reshape(trajectories(iTrajectory).x(ti), [], 1);
                yCell{iTrajectory} = reshape(trajectories(iTrajectory).y(ti), [], 1);
                xDotCell{iTrajectory} = reshape(trajectories(iTrajectory).u(ti), [], 1);
                yDotCell{iTrajectory} = reshape(trajectories(iTrajectory).v(ti), [], 1);
            end

            allT = vertcat(tCell{:});
            allX = vertcat(xCell{:});
            allY = vertcat(yCell{:});
            allXDot = vertcat(xDotCell{:});
            allYDot = vertcat(yDotCell{:});
            mxDotAll = reshape(fit.centerOfMassTrajectory.u(allT), [], 1);
            myDotAll = reshape(fit.centerOfMassTrajectory.v(allT), [], 1);
            uMesoscaleObserved = reshape(fit.uMesoscale(allT, allX, allY), [], 1);
            vMesoscaleObserved = reshape(fit.vMesoscale(allT, allX, allY), [], 1);
            B = reshape(BSpline.matrixForDataPoints(allT, knotPoints=fit.fastKnotPoints, S=fit.fastS), numel(allT), []);
            uMesoscaleComObserved = B * (B \ uMesoscaleObserved);
            vMesoscaleComObserved = B * (B \ vMesoscaleObserved);
            uBackgroundObserved = mxDotAll - uMesoscaleComObserved;
            vBackgroundObserved = myDotAll - vMesoscaleComObserved;
            uSubmesoscaleObserved = allXDot - uBackgroundObserved - uMesoscaleObserved;
            vSubmesoscaleObserved = allYDot - vBackgroundObserved - vMesoscaleObserved;

            kappaValues = zeros(nTrajectories, 1);
            sampleStartIndex = 1;

            for iTrajectory = 1:nTrajectories
                ti = tCell{iTrajectory};
                duration = ti(end) - ti(1);
                nSamples = numel(ti);
                sampleIndices = sampleStartIndex:(sampleStartIndex + nSamples - 1);
                componentS = min(3, nSamples - 1);
                uSpline = InterpolatingSpline.fromGriddedValues(ti, uSubmesoscaleObserved(sampleIndices), S=componentS);
                vSpline = InterpolatingSpline.fromGriddedValues(ti, vSubmesoscaleObserved(sampleIndices), S=componentS);
                xEnd = cumsum(uSpline);
                yEnd = cumsum(vSpline);
                kappaValues(iTrajectory) = (xEnd(ti(end)).^2 + yEnd(ti(end)).^2) / (4 * duration);
                sampleStartIndex = sampleStartIndex + nSamples;
            end

            kappaEstimate = mean(kappaValues);
        end

        function scores = directConsensusScores(summary, scoreIndices, mesoscaleConstraint)
            nBootstraps = size(summary.uCenter, 2);
            scores = struct( ...
                "uv", zeros(1, nBootstraps), ...
                "strain", zeros(1, nBootstraps), ...
                "zeta", zeros(1, nBootstraps), ...
                "joint", zeros(1, nBootstraps));

            for iScore = reshape(scoreIndices, 1, [])
                scores.uv = scores.uv + GriddedStreamfunctionBootstrapUnitTests.directScore2DBlock( ...
                    summary.uCenter(iScore, :), summary.vCenter(iScore, :));

                if mesoscaleConstraint ~= "zeroStrain"
                    scores.strain = scores.strain + GriddedStreamfunctionBootstrapUnitTests.directScore2DBlock( ...
                        summary.sigma_n(iScore, :), summary.sigma_s(iScore, :));
                end

                if mesoscaleConstraint ~= "zeroVorticity"
                    scores.zeta = scores.zeta + GriddedStreamfunctionBootstrapUnitTests.directScore1DBlock( ...
                        summary.zeta(iScore, :));
                end
            end

            scores.joint = scores.uv + scores.strain + scores.zeta;
        end

        function score = directScore2DBlock(xValues, yValues)
            xValues = reshape(xValues, [], 1);
            yValues = reshape(yValues, [], 1);
            score = zeros(1, numel(xValues));
            [isConstantX, isConstantY] = GriddedStreamfunctionBootstrapUnitTests.constantFlags(xValues, yValues);

            if isConstantX && isConstantY
                return
            end

            if isConstantX
                score = GriddedStreamfunctionBootstrapUnitTests.directScore1DBlock(yValues);
                return
            end

            if isConstantY
                score = GriddedStreamfunctionBootstrapUnitTests.directScore1DBlock(xValues);
                return
            end

            model = KernelDensityEstimate.fromData([xValues, yValues]);
            pointDensity = GriddedStreamfunctionBootstrapUnitTests.directGaussianDensity2D(model.data, model.bandwidth, [xValues, yValues]);
            score = log10(GriddedStreamfunctionBootstrapUnitTests.safeDensity(pointDensity)).';
        end

        function score = directScore1DBlock(values)
            values = reshape(values, [], 1);
            score = zeros(1, numel(values));

            if GriddedStreamfunctionBootstrapUnitTests.isConstant(values)
                return
            end

            model = KernelDensityEstimate.fromData(values);
            pointDensity = GriddedStreamfunctionBootstrapUnitTests.directGaussianDensity1D(model.data, model.bandwidth, values);
            score = log10(GriddedStreamfunctionBootstrapUnitTests.safeDensity(pointDensity)).';
        end

        function density = safeDensity(values)
            density = reshape(values, [], 1);
            density(~isfinite(density) | density <= 0) = realmin("double");
        end

        function [isConstantX, isConstantY] = constantFlags(xValues, yValues)
            isConstantX = GriddedStreamfunctionBootstrapUnitTests.isConstant(xValues);
            isConstantY = GriddedStreamfunctionBootstrapUnitTests.isConstant(yValues);
        end

        function tf = isConstant(values)
            values = reshape(values, [], 1);
            scale = max([1; abs(values)]);
            tf = (max(values) - min(values)) <= 1e-12 * scale;
        end

        function density = directGaussianDensity1D(data, bandwidth, queryPoints)
            queryPoints = reshape(queryPoints, [], 1);
            density = mean(exp(-0.5 * ((queryPoints - data.')/bandwidth).^2), 2)/(bandwidth * sqrt(2*pi));
        end

        function density = directGaussianDensity2D(data, bandwidth, queryPoints)
            density = mean(exp(-0.5 * (((queryPoints(:, 1) - data(:, 1).')/bandwidth(1)).^2 + ((queryPoints(:, 2) - data(:, 2).')/bandwidth(2)).^2)), 2) ...
                /(2*pi * bandwidth(1) * bandwidth(2));
        end
    end
end
