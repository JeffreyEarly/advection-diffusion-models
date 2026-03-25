classdef IntegratorUnitTests < matlab.unittest.TestCase

    methods (Test)
        function rk4MatchesExponentialGrowth(testCase)
            dt = 0.01;
            t = (0:dt:1).';

            integrator = Integrator(@(~, y) y, 1, dt);
            y = zeros(size(t));
            y(1) = integrator.currentY;
            for iTime = 2:length(t)
                y(iTime) = integrator.StepForwardOneIncrement;
            end

            testCase.verifyLessThanOrEqual(max(abs(y - exp(t))), 1e-8);
        end

        function diffusivityIntegratorKeepsParticlesInsideBounds(testCase)
            rng(1,"twister");

            nParticles = 512;
            dt = 0.25;
            t = (0:dt:5).';
            lowerBound = -1.5;
            upperBound = 1.5;

            integrator = IntegratorWithDiffusivity( ...
                @(~, y) zeros(size(y)), ...
                zeros(nParticles,1), ...
                dt, ...
                0.75, ...
                lowerBound, ...
                upperBound);

            y = integrator.IntegrateAlongDimension(t);
            y = squeeze(y(:,1,:));

            testCase.verifyGreaterThanOrEqual(min(y,[],"all"), lowerBound - 10*eps(lowerBound));
            testCase.verifyLessThanOrEqual(max(y,[],"all"), upperBound + 10*eps(upperBound));
        end

        function diffusivityIntegratorMatchesExpectedKappa(testCase)
            rng(1,"twister");

            nParticles = 20000;
            nDims = 2;
            dt = 0.25;
            nSteps = 160;
            totalTime = nSteps*dt;
            kappa = 0.75;

            integrator = IntegratorWithDiffusivity( ...
                @(~, y) zeros(size(y)), ...
                zeros(nParticles,nDims), ...
                dt, ...
                kappa, ...
                [-Inf -Inf], ...
                [Inf Inf]);

            y = integrator.currentY;
            for iStep = 1:nSteps
                y = integrator.StepForwardOneIncrement;
            end

            dimensionVariance = var(y,1,1);
            kappaEstimate = mean(sum(y.^2,2))/(2*nDims*totalTime);

            testCase.verifyEqual(kappaEstimate, kappa, "RelTol", 0.05);
            testCase.verifyEqual(dimensionVariance, 2*kappa*totalTime*ones(1,nDims), "RelTol", 0.05);
            testCase.verifyEqual(mean(y,1), zeros(1,nDims), "AbsTol", 0.2);
        end

        function diffusivityIntegratorMatchesReflectedWallStatistics(testCase)
            rng(1,"twister");

            nParticles = 20000;
            dt = 0.25;
            kappa = 0.75;

            integrator = IntegratorWithDiffusivity( ...
                @(~, y) zeros(size(y)), ...
                zeros(nParticles,1), ...
                dt, ...
                kappa, ...
                0, ...
                Inf);

            y = integrator.StepForwardOneIncrement;

            expectedMean = 2*sqrt(kappa*dt/pi);
            expectedSecondMoment = 2*kappa*dt;
            expectedVariance = expectedSecondMoment*(1 - 2/pi);

            testCase.verifyGreaterThanOrEqual(min(y), 0);
            testCase.verifyEqual(mean(y), expectedMean, "RelTol", 0.05);
            testCase.verifyEqual(mean(y.^2), expectedSecondMoment, "RelTol", 0.05);
            testCase.verifyEqual(var(y,1), expectedVariance, "RelTol", 0.05);
        end

        function eulerMaruyamaMatchesExpectedKappa(testCase)
            rng(1,"twister");

            nParticles = 20000;
            nDims = 2;
            kappa = 0.75;
            dt = 0.25;
            totalTime = 40;
            t = (0:dt:totalTime).';

            integrator = IntegratorEulerMaruyama( ...
                @(~, y) zeros(size(y)), ...
                @(~, y) sqrt(2*kappa)*ones(size(y)), ...
                zeros(nParticles,nDims), ...
                dt);

            y = integrator.IntegrateAlongDimension(t);
            y = squeeze(y(:,:,end));

            dimensionVariance = var(y,1,1);
            kappaEstimate = mean(sum(y.^2,2))/(2*nDims*totalTime);

            testCase.verifyEqual(kappaEstimate, kappa, "RelTol", 0.05);
            testCase.verifyEqual(dimensionVariance, 2*kappa*totalTime*ones(1,nDims), "RelTol", 0.05);
            testCase.verifyEqual(mean(y,1), zeros(1,nDims), "AbsTol", 0.2);
        end
    end

end
