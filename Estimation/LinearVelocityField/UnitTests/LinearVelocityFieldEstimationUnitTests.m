classdef LinearVelocityFieldEstimationUnitTests < matlab.unittest.TestCase

    methods (Test)
        function estimateConstantLinearVelocityFieldParameters(testCase)
            sigma = 4.0e-6;
            theta = pi/9;
            zeta = -1.5e-6;
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
            integrator = AdvectionDiffusionIntegrator(model,0);

            [x0, y0] = ndgrid([-600 0 600],[-400 0 400]);
            [t, x, y] = integrator.particleTrajectories(x0(:),y0(:),12*3600,900);

            parameters = EstimateLinearVelocityFieldParameters(x, y, t, [ModelParameter.u0v0, ModelParameter.strain, ModelParameter.vorticity], 1);

            testCase.verifyEqual(parameters.u0, u0, "AbsTol", 5e-6);
            testCase.verifyEqual(parameters.v0, v0, "AbsTol", 5e-6);
            testCase.verifyEqual(parameters.sigma_n, model.sigma_n, "AbsTol", 5e-10);
            testCase.verifyEqual(parameters.sigma_s, model.sigma_s, "AbsTol", 5e-10);
            testCase.verifyEqual(parameters.zeta, zeta, "AbsTol", 5e-10);
        end

        function estimateTimeVaryingBackgroundVelocityWithSplines(testCase)
            t = linspace(0,12*3600,49).';
            T = t(end);

            u0 = 0.05 + 0.03*(t/T);
            v0 = -0.04 + 0.02*(t/T);

            xCenter = 0.05*t + 0.03*(t.^2)/(2*T);
            yCenter = -0.04*t + 0.02*(t.^2)/(2*T);

            [xOffset, yOffset] = ndgrid([-300 0 300],[-200 0 200]);
            x = xCenter + reshape(xOffset,1,[]);
            y = yCenter + reshape(yOffset,1,[]);

            parameters = EstimateLinearVelocityFieldParameters(x, y, t, ModelParameter.u0v0, 4);

            testCase.verifySize(parameters.u0,size(t));
            testCase.verifySize(parameters.v0,size(t));
            testCase.verifyLessThanOrEqual(max(abs(parameters.u0 - u0)), 5e-7);
            testCase.verifyLessThanOrEqual(max(abs(parameters.v0 - v0)), 5e-7);
            testCase.verifyEqual(parameters.sigma_n, zeros(size(t)), "AbsTol", 1e-12);
            testCase.verifyEqual(parameters.sigma_s, zeros(size(t)), "AbsTol", 1e-12);
            testCase.verifyEqual(parameters.zeta, zeros(size(t)), "AbsTol", 1e-12);
        end
    end

end
