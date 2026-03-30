classdef TensorSplineStreamfunctionUnitTests < matlab.unittest.TestCase

    methods (Test)
        function synchronousLinearFieldRecovery(testCase)
            [model, t, x, y] = TensorSplineStreamfunctionUnitTests.synchronousLinearFieldData();

            fit = TensorSplineStreamfunction.fitFromSynchronousTrajectories(x, y, t);

            xcTrue = mean(x, 2);
            ycTrue = mean(y, 2);
            tGrid = repmat(t, 1, size(x, 2));
            uExpected = model.u(t, xcTrue, ycTrue);
            vExpected = model.v(t, xcTrue, ycTrue);

            testCase.verifyLessThanOrEqual(max(abs(fit.xCenter(t) - xcTrue)), 2e-1);
            testCase.verifyLessThanOrEqual(max(abs(fit.yCenter(t) - ycTrue)), 2e-1);
            testCase.verifyLessThanOrEqual(max(abs(fit.uBackground(t) - uExpected)), 2e-4);
            testCase.verifyLessThanOrEqual(max(abs(fit.vBackground(t) - vExpected)), 2e-4);
            testCase.verifyLessThanOrEqual(max(abs(fit.sigma_n(tGrid, x, y) - model.sigma_n), [], "all"), 2e-8);
            testCase.verifyLessThanOrEqual(max(abs(fit.sigma_s(tGrid, x, y) - model.sigma_s), [], "all"), 2e-8);
            testCase.verifyLessThanOrEqual(max(abs(fit.zeta(tGrid, x, y) - model.zeta), [], "all"), 2e-8);
            testCase.verifyLessThanOrEqual(max(abs(fit.u(tGrid, x, y) - model.u(tGrid, x, y)), [], "all"), 2e-4);
            testCase.verifyLessThanOrEqual(max(abs(fit.v(tGrid, x, y) - model.v(tGrid, x, y)), [], "all"), 2e-4);
        end

        function asynchronousLinearFieldRecovery(testCase)
            sigma = 4.0e-6;
            theta = pi/9;
            zeta = -1.5e-6;
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
            integrator = AdvectionDiffusionIntegrator(model, 0);

            [x0, y0] = ndgrid([-600 0 600], [-400 0 400]);
            [tFine, xFine, yFine] = integrator.particleTrajectories(x0(:), y0(:), 12*3600, 300);

            xcFine = mean(xFine, 2);
            ycFine = mean(yFine, 2);
            [xCell, yCell, tCell, xcCell, ycCell] = TensorSplineStreamfunctionUnitTests.asynchronousCells(tFine, xFine, yFine, xcFine, ycFine);

            for iDrifter = 1:numel(tCell)
                testCase.verifyGreaterThan(numel(unique(diff(tCell{iDrifter}))), 1);
            end

            fit = TensorSplineStreamfunction.fitFromTrajectories(xCell, yCell, tCell);

            uBackgroundError = zeros(numel(tCell), 1);
            vBackgroundError = zeros(numel(tCell), 1);
            uTotalError = zeros(numel(tCell), 1);
            vTotalError = zeros(numel(tCell), 1);
            sigmaNError = zeros(numel(tCell), 1);
            sigmaSError = zeros(numel(tCell), 1);
            zetaError = zeros(numel(tCell), 1);

            for iDrifter = 1:numel(tCell)
                ti = tCell{iDrifter};
                xi = xCell{iDrifter};
                yi = yCell{iDrifter};
                xci = xcCell{iDrifter};
                yci = ycCell{iDrifter};

                uBackgroundError(iDrifter) = max(abs(fit.uBackground(ti) - model.u(ti, xci, yci)));
                vBackgroundError(iDrifter) = max(abs(fit.vBackground(ti) - model.v(ti, xci, yci)));
                uTotalError(iDrifter) = max(abs(fit.u(ti, xi, yi) - model.u(ti, xi, yi)));
                vTotalError(iDrifter) = max(abs(fit.v(ti, xi, yi) - model.v(ti, xi, yi)));
                sigmaNError(iDrifter) = max(abs(fit.sigma_n(ti, xi, yi) - model.sigma_n));
                sigmaSError(iDrifter) = max(abs(fit.sigma_s(ti, xi, yi) - model.sigma_s));
                zetaError(iDrifter) = max(abs(fit.zeta(ti, xi, yi) - model.zeta));
            end

            testCase.verifyLessThanOrEqual(max(uBackgroundError), 3e-4);
            testCase.verifyLessThanOrEqual(max(vBackgroundError), 3e-4);
            testCase.verifyLessThanOrEqual(max(uTotalError), 3e-4);
            testCase.verifyLessThanOrEqual(max(vTotalError), 3e-4);
            testCase.verifyLessThanOrEqual(max(sigmaNError), 3e-8);
            testCase.verifyLessThanOrEqual(max(sigmaSError), 3e-8);
            testCase.verifyLessThanOrEqual(max(zetaError), 3e-8);
        end

        function synchronousWrapperMatchesCellFit(testCase)
            [~, t, x, y] = TensorSplineStreamfunctionUnitTests.synchronousLinearFieldData();

            fitMatrix = TensorSplineStreamfunction.fitFromSynchronousTrajectories(x, y, t);

            xCell = mat2cell(x, size(x, 1), ones(1, size(x, 2)));
            yCell = mat2cell(y, size(y, 1), ones(1, size(y, 2)));
            tCell = repmat({t}, size(x, 2), 1);
            fitCell = TensorSplineStreamfunction.fitFromTrajectories(xCell, yCell, tCell);

            testCase.verifyEqual(fitMatrix.centerXSpline.xi, fitCell.centerXSpline.xi, "AbsTol", 1e-12);
            testCase.verifyEqual(fitMatrix.centerYSpline.xi, fitCell.centerYSpline.xi, "AbsTol", 1e-12);
            testCase.verifyEqual(fitMatrix.streamfunctionSpline.xi, fitCell.streamfunctionSpline.xi, "AbsTol", 1e-12);
        end

        function explicitPsiKnotPointsMatchDefaultConstruction(testCase)
            [~, t, x, y] = TensorSplineStreamfunctionUnitTests.synchronousLinearFieldData();

            fitDefault = TensorSplineStreamfunction.fitFromSynchronousTrajectories(x, y, t);
            fitExplicit = TensorSplineStreamfunction.fitFromSynchronousTrajectories(x, y, t, ...
                psiS=fitDefault.streamfunctionSpline.S, psiKnotPoints=fitDefault.streamfunctionSpline.knotPoints);

            testCase.verifyEqual(fitDefault.centerXSpline.xi, fitExplicit.centerXSpline.xi, "AbsTol", 1e-12);
            testCase.verifyEqual(fitDefault.centerYSpline.xi, fitExplicit.centerYSpline.xi, "AbsTol", 1e-12);
            testCase.verifyEqual(fitDefault.streamfunctionSpline.xi, fitExplicit.streamfunctionSpline.xi, "AbsTol", 1e-12);
        end

        function evaluationInputShapes(testCase)
            [~, t, x, y] = TensorSplineStreamfunctionUnitTests.synchronousLinearFieldData();
            fit = TensorSplineStreamfunction.fitFromSynchronousTrajectories(x, y, t);

            tGrid = repmat(t, 1, size(x, 2));
            scalarTime = t(4);
            scalarRow = repmat(scalarTime, 1, size(x, 2));

            testCase.verifyEqual(fit.u(scalarTime, x(4, :), y(4, :)), fit.u(scalarRow, x(4, :), y(4, :)), "AbsTol", 1e-12);
            testCase.verifyEqual(fit.psi(t, x, y), fit.psi(tGrid, x, y), "AbsTol", 1e-12);
            testCase.verifyEqual(fit.u(t, x, y), fit.u(tGrid, x, y), "AbsTol", 1e-12);
            testCase.verifyEqual(fit.v(t, x, y), fit.v(tGrid, x, y), "AbsTol", 1e-12);

            testCase.verifyError(@() fit.u(t, x(1:end-1, :), y), ...
                "TensorSplineStreamfunction:EvaluationSizeMismatch");
            testCase.verifyError(@() fit.u(t(1:end-1), x, y), ...
                "TensorSplineStreamfunction:EvaluationTimeSizeMismatch");
        end
    end

    methods (Static, Access = private)
        function [model, t, x, y] = synchronousLinearFieldData()
            sigma = 4.0e-6;
            theta = pi/9;
            zeta = -1.5e-6;
            u0 = 0.08;
            v0 = -0.05;

            model = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
            integrator = AdvectionDiffusionIntegrator(model, 0);

            [x0, y0] = ndgrid([-600 0 600], [-400 0 400]);
            [t, x, y] = integrator.particleTrajectories(x0(:), y0(:), 12*3600, 900);
        end

        function [xCell, yCell, tCell, xcCell, ycCell] = asynchronousCells(tFine, xFine, yFine, xcFine, ycFine)
            nDrifters = size(xFine, 2);
            xCell = cell(nDrifters, 1);
            yCell = cell(nDrifters, 1);
            tCell = cell(nDrifters, 1);
            xcCell = cell(nDrifters, 1);
            ycCell = cell(nDrifters, 1);

            % Opposite drifters share the same offset so each pseudo-center
            % sample remains balanced while individual drifters still have
            % nonuniform time intervals.
            drifterOffset = [-1; 0; 1; 1; 0; 1; 1; 0; -1];
            baseIndices = 3:3:(numel(tFine) - 2);
            jitterPattern = [0; 1; 0; -1; 0; 1; 0; -1; 0; 0];
            jitter = repmat(jitterPattern, ceil(numel(baseIndices) / numel(jitterPattern)), 1);
            jitter = jitter(1:numel(baseIndices));

            for iDrifter = 1:nDrifters
                indices = baseIndices + transpose(jitter) + drifterOffset(iDrifter);

                xCell{iDrifter} = xFine(indices, iDrifter);
                yCell{iDrifter} = yFine(indices, iDrifter);
                tCell{iDrifter} = tFine(indices);
                xcCell{iDrifter} = xcFine(indices);
                ycCell{iDrifter} = ycFine(indices);
            end
        end
    end
end
