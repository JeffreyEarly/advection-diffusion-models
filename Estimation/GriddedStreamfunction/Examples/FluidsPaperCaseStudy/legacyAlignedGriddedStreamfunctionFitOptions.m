function options = legacyAlignedGriddedStreamfunctionFitOptions(t, x, y, dof)
arguments
    t (:,1) double {mustBeReal, mustBeFinite}
    x (:,:) double {mustBeReal, mustBeFinite}
    y (:,:) double {mustBeReal, mustBeFinite}
    dof (1,1) double {mustBeInteger, mustBePositive} = 4
end

if ~isequal(size(x), size(y))
    error("GriddedStreamfunction:ExampleCoordinateSizeMismatch", ...
        "x and y must have the same size.");
end
if size(x, 1) ~= numel(t)
    error("GriddedStreamfunction:ExampleTimeSizeMismatch", ...
        "The first dimension of x and y must match numel(t).");
end

[timeBasisMatrix, tKnotPoints, timeSplineDegree] = legacyTimeBasis(t, dof);
psiS = [2 2 timeSplineDegree];

mx = timeBasisMatrix * (timeBasisMatrix \ mean(x, 2));
my = timeBasisMatrix * (timeBasisMatrix \ mean(y, 2));
q = x - mx;
r = y - my;

options = struct( ...
    "dof", dof, ...
    "timeBasisMatrix", timeBasisMatrix, ...
    "timeSplineDegree", timeSplineDegree, ...
    "fastKnotPoints", tKnotPoints, ...
    "fastS", timeSplineDegree, ...
    "psiKnotPoints", {{ ...
        paddedQuadraticDomain(q) ...
        paddedQuadraticDomain(r) ...
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
