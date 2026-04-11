classdef KernelDensityEstimateUnitTests < matlab.unittest.TestCase

    methods (Test)
        function fitOneDimensionalDataMatchesLegacyBandwidth(testCase)
            data = KernelDensityEstimateUnitTests.oneDimensionalData();
            model = KernelDensityEstimate.fromData(data);

            testCase.verifyClass(model, "KernelDensityEstimate")
            testCase.verifyEqual(model.numDimensions, 1)
            testCase.verifyEqual(model.bandwidth, 2.21333986156179, "AbsTol", 1e-12)
        end

        function fitOneDimensionalDataMatchesFrozenBandwidthOnAdditionalDataset(testCase)
            datasets = KernelDensityEstimateUnitTests.oneDimensionalDatasets();
            model = KernelDensityEstimate.fromData(datasets(2).data);

            testCase.verifyEqual(model.bandwidth, datasets(2).bandwidth, "AbsTol", 1e-12)
        end

        function fitOneDimensionalDataIsPermutationInvariant(testCase)
            datasets = KernelDensityEstimateUnitTests.oneDimensionalDatasets();

            for iDataset = 1:numel(datasets)
                data = datasets(iDataset).data;
                baseline = KernelDensityEstimate.fromData(data);
                permuted = KernelDensityEstimate.fromData(data(end:-1:1));

                testCase.verifyEqual(permuted.bandwidth, baseline.bandwidth, "AbsTol", 1e-12)
            end
        end

        function fitOneDimensionalDataIsTranslationInvariant(testCase)
            datasets = KernelDensityEstimateUnitTests.oneDimensionalDatasets();
            translation = 8;

            for iDataset = 1:numel(datasets)
                baseline = KernelDensityEstimate.fromData(datasets(iDataset).data);
                translated = KernelDensityEstimate.fromData(datasets(iDataset).data + translation);

                testCase.verifyEqual(translated.bandwidth, baseline.bandwidth, "AbsTol", 1e-12)
            end
        end

        function fitOneDimensionalDataIsPositiveScaleEquivariant(testCase)
            datasets = KernelDensityEstimateUnitTests.oneDimensionalDatasets();
            scaleFactor = 4;

            for iDataset = 1:numel(datasets)
                baseline = KernelDensityEstimate.fromData(datasets(iDataset).data);
                scaled = KernelDensityEstimate.fromData(datasets(iDataset).data * scaleFactor);

                testCase.verifyEqual(scaled.bandwidth, baseline.bandwidth * scaleFactor, "AbsTol", 1e-10, "RelTol", 1e-10)
            end
        end

        function fitTwoDimensionalDataPreservesLegacyKde2dBandwidthSmokeTest(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            model = KernelDensityEstimate.fromData(data);

            testCase.verifyEqual(model.numDimensions, 2)
            testCase.verifyEqual(model.bandwidth, [1.78508660817316 1.09168358178069], "AbsTol", 1e-12)
        end

        function fitTwoDimensionalDataMatchesReferenceBandwidthOnFrozenDatasets(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();

            for iDataset = 1:numel(datasets)
                model = KernelDensityEstimate.fromData(datasets(iDataset).data);
                diagnostics = referenceBandwidthForTwoDimensionalData(datasets(iDataset).data);

                testCase.verifyEqual(model.bandwidth, diagnostics.bandwidth, "AbsTol", 1e-10, "RelTol", 1e-10)
            end
        end

        function referenceSolverResidualIsSmallOnFrozenDatasets(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();

            for iDataset = 1:numel(datasets)
                diagnostics = referenceBandwidthForTwoDimensionalData(datasets(iDataset).data);

                testCase.verifyLessThanOrEqual(abs(diagnostics.objectiveResidual), 1e-10)
            end
        end

        function fitTwoDimensionalDataIsPermutationInvariant(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();

            for iDataset = 1:numel(datasets)
                data = datasets(iDataset).data;
                baseline = KernelDensityEstimate.fromData(data);
                permuted = KernelDensityEstimate.fromData(data(end:-1:1, :));

                testCase.verifyEqual(permuted.bandwidth, baseline.bandwidth, "AbsTol", 1e-12)
            end
        end

        function fitTwoDimensionalDataIsTranslationInvariant(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();
            translation = [8 -4];

            for iDataset = 1:numel(datasets)
                baseline = KernelDensityEstimate.fromData(datasets(iDataset).data);
                translated = KernelDensityEstimate.fromData(datasets(iDataset).data + translation);

                testCase.verifyEqual(translated.bandwidth, baseline.bandwidth, "AbsTol", 1e-12)
            end
        end

        function fitTwoDimensionalDataIsComponentwiseScaleInvariant(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();
            scaleFactors = [4 0.25];

            for iDataset = 1:numel(datasets)
                baseline = KernelDensityEstimate.fromData(datasets(iDataset).data);
                scaled = KernelDensityEstimate.fromData(datasets(iDataset).data .* scaleFactors);

                testCase.verifyEqual(scaled.bandwidth, baseline.bandwidth .* scaleFactors, "AbsTol", 1e-10, "RelTol", 1e-10)
            end
        end

        function densityAtMatchesDirectGaussianReferenceInOneDimension(testCase)
            data = KernelDensityEstimateUnitTests.oneDimensionalData();
            model = KernelDensityEstimate.fromData(data);
            queryPoints = linspace(-4, 4, 11).';

            expected = KernelDensityEstimateUnitTests.directGaussianDensity1D(model.data, model.bandwidth, queryPoints);

            testCase.verifyEqual(model.densityAt(queryPoints), expected, "AbsTol", 1e-12)
            testCase.verifyEqual(model.densityAt(queryPoints.'), expected.', "AbsTol", 1e-12)
        end

        function densityAtMatchesDirectGaussianReferenceInTwoDimensions(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            model = KernelDensityEstimate.fromData(data);
            queryPoints = [-1.5 -0.8; -0.25 0.5; 1.6 1.2];

            expected = KernelDensityEstimateUnitTests.directGaussianDensity2D(model.data, model.bandwidth, queryPoints);

            testCase.verifyEqual(model.densityAt(queryPoints), expected, "AbsTol", 1e-12)
        end

        function densityOnGridMatchesPointEvaluationAndSupportsArbitraryGridSizes(testCase)
            model1D = KernelDensityEstimate.fromData(KernelDensityEstimateUnitTests.oneDimensionalData());
            [density1D, gridVectors1D] = model1D.densityOnGrid(gridSize=300);

            testCase.verifySize(density1D, [300 1])
            testCase.verifyEqual(gridVectors1D{1}(1), model1D.minimum, "AbsTol", 1e-12)
            testCase.verifyEqual(gridVectors1D{1}(end), model1D.maximum, "AbsTol", 1e-12)
            testCase.verifyEqual(density1D, model1D.densityAt(gridVectors1D{1}), "AbsTol", 1e-12)

            model2D = KernelDensityEstimate.fromData(KernelDensityEstimateUnitTests.twoDimensionalData());
            [density2D, gridVectors2D] = model2D.densityOnGrid(gridSize=[192 320]);
            [X, Y] = ndgrid(gridVectors2D{1}, gridVectors2D{2});

            testCase.verifySize(density2D, [192 320])
            testCase.verifyEqual(density2D, reshape(model2D.densityAt([X(:), Y(:)]), size(X)), "AbsTol", 1e-12)
        end

        function cdfOnGridMatchesDirectGaussianReference(testCase)
            data = KernelDensityEstimateUnitTests.oneDimensionalData();
            model = KernelDensityEstimate.fromData(data);
            [cdfValues, gridVectors] = model.cdfOnGrid(gridSize=301);

            expected = KernelDensityEstimateUnitTests.directGaussianCdf1D(model.data, model.bandwidth, gridVectors{1});

            testCase.verifyEqual(cdfValues, expected, "AbsTol", 1e-12)
            testCase.verifyGreaterThanOrEqual(min(cdfValues), 0)
            testCase.verifyLessThanOrEqual(max(cdfValues), 1)
        end

        function bandwidthGridSizeRoundsUpInternally(testCase)
            data = KernelDensityEstimateUnitTests.oneDimensionalData();
            roundedModel = KernelDensityEstimate.fromData(data, bandwidthGridSize=300);
            exactModel = KernelDensityEstimate.fromData(data, bandwidthGridSize=512);

            testCase.verifyEqual(roundedModel.bandwidth, exactModel.bandwidth, "AbsTol", 1e-12)
        end

        function twoDimensionalBandwidthGridSizeRoundsUpInternally(testCase)
            datasets = KernelDensityEstimateUnitTests.twoDimensionalDatasets();

            for iDataset = 1:numel(datasets)
                roundedModel = KernelDensityEstimate.fromData(datasets(iDataset).data, bandwidthGridSize=300);
                exactModel = KernelDensityEstimate.fromData(datasets(iDataset).data, bandwidthGridSize=512);

                testCase.verifyEqual(roundedModel.bandwidth, exactModel.bandwidth, "AbsTol", 1e-12)
            end
        end

        function invalidDimensionsAndBoundsRaiseStructuredErrors(testCase)
            data = KernelDensityEstimateUnitTests.oneDimensionalData();

            testCase.verifyError(@() KernelDensityEstimate.fromData(rand(5, 3)), "KernelDensityEstimate:UnsupportedDimension")
            testCase.verifyError(@() KernelDensityEstimate.fromData(data, minimum=1, maximum=0), "KernelDensityEstimate:InvalidBounds")
            testCase.verifyError(@() KernelDensityEstimate.fromData(data, minimum=[-1 -1], maximum=[1 1]), "KernelDensityEstimate:InvalidBounds")
        end

        function densityLevelForCdfMatchesLegacyContourHelper(testCase)
            data = KernelDensityEstimateUnitTests.twoDimensionalData();
            model = KernelDensityEstimate.fromData(data);
            [density, gridVectors] = model.densityOnGrid(gridSize=[64 64]);
            pctTarget = [0.1; 0.5; 0.8];

            expected = KernelDensityEstimateUnitTests.legacyDensityLevelForCDF(gridVectors, density, pctTarget);
            actual = DensityLevelForCDF(gridVectors, density, pctTarget);

            testCase.verifyEqual(actual, expected, "AbsTol", 1e-12)
        end
    end

    methods (Static, Access = private)
        function data = oneDimensionalData()
            datasets = KernelDensityEstimateUnitTests.oneDimensionalDatasets();
            data = datasets(1).data;
        end

        function datasets = oneDimensionalDatasets()
            datasets = [ ...
                struct( ...
                    "name", "Original sparse sample", ...
                    "data", [-3.1; -1.2; -0.8; 0.1; 0.4; 0.9; 1.7; 2.4; 3.0], ...
                    "bandwidth", 2.21333986156179), ...
                struct( ...
                    "name", "Bimodal cloud with repeated sample", ...
                    "data", [-5.0; -4.5; -4.0; -3.0; 2.5; 3.0; 3.0; 4.0; 4.5; 5.0], ...
                    "bandwidth", 1.73785891998607)];
        end

        function data = twoDimensionalData()
            data = [-2.0 -1.1; -1.4 -0.2; -0.8 0.4; -0.2 0.9; 0.5 -0.6; 1.0 0.2; 1.7 1.1; 2.4 1.8];
        end

        function datasets = twoDimensionalDatasets()
            datasets = [ ...
                struct( ...
                    "name", "Compact anisotropic cloud", ...
                    "data", [-2.00 -1.00; -1.50 -0.25; -1.00 0.50; -0.25 1.00; 0.50 -0.50; 1.00 0.25; 1.75 1.00; 2.50 1.75]), ...
                struct( ...
                    "name", "Well-separated axis-aligned bimodal cloud", ...
                    "data", [-5.0 -1.50; -4.5 -1.00; -4.0 -1.25; -3.5 -0.50; -3.0 -1.00; 2.5 1.50; 3.0 2.00; 3.5 1.75; 4.0 2.25; 4.5 2.00]), ...
                struct( ...
                    "name", "Strongly anisotropic single cloud", ...
                    "data", [-6.00 0.12500; -4.75 -0.12500; -3.50 0.06250; -2.25 -0.06250; -1.00 0.03125; 0.25 -0.03125; 1.50 0.09375; 2.75 -0.09375; 4.00 0.046875; 5.50 -0.046875])];
        end

        function density = directGaussianDensity1D(data, bandwidth, queryPoints)
            queryPoints = reshape(queryPoints, [], 1);
            density = mean(exp(-0.5 * ((queryPoints - data.')/bandwidth).^2), 2)/(bandwidth * sqrt(2*pi));
        end

        function density = directGaussianDensity2D(data, bandwidth, queryPoints)
            density = mean(exp(-0.5 * (((queryPoints(:, 1) - data(:, 1).')/bandwidth(1)).^2 + ((queryPoints(:, 2) - data(:, 2).')/bandwidth(2)).^2)), 2) ...
                /(2*pi * bandwidth(1) * bandwidth(2));
        end

        function cdfValues = directGaussianCdf1D(data, bandwidth, queryPoints)
            queryPoints = reshape(queryPoints, [], 1);
            cdfValues = mean(0.5 * erfc(-(queryPoints - data.')/(bandwidth * sqrt(2))), 2);
        end

        function dLevels = legacyDensityLevelForCDF(gridVectors, density, pctTarget)
            [X, Y] = meshgrid(gridVectors{1}, gridVectors{2});
            densityMesh = density.';
            nLevels = 25;
            level = zeros(nLevels, 1);
            pctEnclosed = zeros(nLevels, 1);
            sigmaNAxis = X(1, :).';
            sigmaSAxis = Y(:, 1);
            contours = contourc(sigmaNAxis, sigmaSAxis, densityMesh, nLevels);
            iContour = 1;
            iLevel = 1;

            while iContour < size(contours, 2)
                level(iLevel) = contours(1, iContour);
                nVertices = contours(2, iContour);
                polygon.Vertices = [contours(1, (iContour + 1):(iContour + nVertices - 1)).', ...
                    contours(2, (iContour + 1):(iContour + nVertices - 1)).'];
                polygon.dx = polygon.Vertices(2:end, 1) - polygon.Vertices(1:end - 1, 1);
                polygon.dy = polygon.Vertices(2:end, 2) - polygon.Vertices(1:end - 1, 2);
                mask = reshape(KernelDensityEstimateUnitTests.isInteriorLegacy(polygon, X(:), Y(:)), size(X));
                pctEnclosed(iLevel) = trapz(sigmaSAxis, trapz(sigmaNAxis, densityMesh .* mask));
                iContour = iContour + nVertices + 1;
                iLevel = iLevel + 1;
            end

            dLevels = interp1(pctEnclosed(1:iLevel - 1), level(1:iLevel - 1), pctTarget);
        end

        function inside = isInteriorLegacy(polygon, x, y)
            windingNumber = zeros(size(x));
            for iVertex = 1:(length(polygon.Vertices) - 1)
                isBelow = polygon.Vertices(iVertex, 2) <= y;
                upwardCrossing = isBelow & polygon.Vertices(iVertex + 1, 2) > y;
                if any(upwardCrossing)
                    windingNumber(upwardCrossing) = windingNumber(upwardCrossing) + ...
                        (KernelDensityEstimateUnitTests.isLeftLegacy(polygon, iVertex, x(upwardCrossing), y(upwardCrossing)) > 0);
                end

                downwardCrossing = ~isBelow & polygon.Vertices(iVertex + 1, 2) <= y;
                if any(downwardCrossing)
                    windingNumber(downwardCrossing) = windingNumber(downwardCrossing) - ...
                        (KernelDensityEstimateUnitTests.isLeftLegacy(polygon, iVertex, x(downwardCrossing), y(downwardCrossing)) < 0);
                end
            end

            inside = abs(windingNumber) > 0;
        end

        function isLeftValue = isLeftLegacy(polygon, iVertex, x, y)
            isLeftValue = polygon.dx(iVertex) .* (y - polygon.Vertices(iVertex, 2)) - (x - polygon.Vertices(iVertex, 1)) .* polygon.dy(iVertex);
        end
    end
end
