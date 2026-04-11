function dLevels = DensityLevelForCDF(gridVectors, density, pctTarget)
% DensityLevelForCDF Return contour levels that enclose target density mass.
%
% `gridVectors` should be the cell-array grid returned by
% `KernelDensityEstimate.densityOnGrid(...)`, with `density` stored in the
% same `ndgrid` ordering. `pctTarget` is a vector of target enclosed-mass
% fractions on `[0, 1]`.

arguments
    gridVectors (1,2) cell
    density (:,:) double {mustBeReal, mustBeFinite, mustBeNonempty}
    pctTarget (:,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(pctTarget, 0), mustBeLessThanOrEqual(pctTarget, 1)}
end

xVector = reshape(gridVectors{1}, [], 1);
yVector = reshape(gridVectors{2}, [], 1);
if ~isnumeric(xVector) || ~isnumeric(yVector) || any(~isfinite(xVector)) || any(~isfinite(yVector))
    error("DensityLevelForCDF:InvalidGridVectors", ...
        "gridVectors must contain finite numeric coordinate vectors.");
end
if ~isequal(size(density), [numel(xVector) numel(yVector)])
    error("DensityLevelForCDF:SizeMismatch", ...
        "density must have size [numel(gridVectors{1}) numel(gridVectors{2})].");
end

nLevels = 25;
level = zeros(nLevels, 1);
pctEnclosed = zeros(nLevels, 1);
[X, Y] = ndgrid(xVector, yVector);
contours = contourc(xVector, yVector, density.', nLevels);
iContour = 1;
iLevel = 1;

while iContour < size(contours, 2)
    level(iLevel) = contours(1, iContour);
    nVertices = contours(2, iContour);
    polygon.Vertices = [contours(1, (iContour + 1):(iContour + nVertices - 1)).', ...
        contours(2, (iContour + 1):(iContour + nVertices - 1)).'];
    polygon.dx = polygon.Vertices(2:end, 1) - polygon.Vertices(1:end - 1, 1);
    polygon.dy = polygon.Vertices(2:end, 2) - polygon.Vertices(1:end - 1, 2);
    mask = reshape(isInterior(polygon, X(:), Y(:)), size(X));
    pctEnclosed(iLevel) = trapz(yVector, trapz(xVector, density .* mask, 1), 2);

    iContour = iContour + nVertices + 1;
    iLevel = iLevel + 1;
end

dLevels = interp1(pctEnclosed(1:iLevel - 1), level(1:iLevel - 1), pctTarget);
end

function isLeftValue = isLeft(polygon, iVertex, x, y)
isLeftValue = polygon.dx(iVertex) .* (y - polygon.Vertices(iVertex, 2)) - (x - polygon.Vertices(iVertex, 1)) .* polygon.dy(iVertex);
end

function inside = isInterior(polygon, x, y)
windingNumber = zeros(size(x));
for iVertex = 1:(length(polygon.Vertices) - 1)
    isBelow = polygon.Vertices(iVertex, 2) <= y;
    upwardCrossing = isBelow & polygon.Vertices(iVertex + 1, 2) > y;
    if any(upwardCrossing)
        windingNumber(upwardCrossing) = windingNumber(upwardCrossing) + (isLeft(polygon, iVertex, x(upwardCrossing), y(upwardCrossing)) > 0);
    end

    downwardCrossing = ~isBelow & polygon.Vertices(iVertex + 1, 2) <= y;
    if any(downwardCrossing)
        windingNumber(downwardCrossing) = windingNumber(downwardCrossing) - (isLeft(polygon, iVertex, x(downwardCrossing), y(downwardCrossing)) < 0);
    end
end

inside = abs(windingNumber) > 0;
end
