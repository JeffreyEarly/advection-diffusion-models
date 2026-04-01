function fitTrajectorySplines(self, trajectories, psiKnotPoints, psiS, fastKnotPoints, fastS)
nTrajectories = numel(trajectories);
tCell = cell(nTrajectories, 1);
xCell = cell(nTrajectories, 1);
yCell = cell(nTrajectories, 1);
xDotCell = cell(nTrajectories, 1);
yDotCell = cell(nTrajectories, 1);

for iTrajectory = 1:nTrajectories
    ti = reshape(trajectories(iTrajectory).t, [], 1);
    xi = reshape(trajectories(iTrajectory).x(ti), [], 1);
    yi = reshape(trajectories(iTrajectory).y(ti), [], 1);
    xDoti = reshape(trajectories(iTrajectory).x.valueAtPoints(ti, D=1), [], 1);
    yDoti = reshape(trajectories(iTrajectory).y.valueAtPoints(ti, D=1), [], 1);

    tCell{iTrajectory} = ti;
    xCell{iTrajectory} = xi;
    yCell{iTrajectory} = yi;
    xDotCell{iTrajectory} = xDoti;
    yDotCell{iTrajectory} = yDoti;
end

allT = vertcat(tCell{:});
allX = vertcat(xCell{:});
allY = vertcat(yCell{:});
allXDot = vertcat(xDotCell{:});
allYDot = vertcat(yDotCell{:});
sampleCounts = reshape(cellfun(@numel, tCell), [], 1);
fitSupportTimes = unique(allT, "sorted");

if isempty(fastKnotPoints)
    representativeTimes = GriddedStreamfunction.representativeObservationTimes(tCell);
    fastSupportTimes = unique(representativeTimes, "sorted");
    if numel(fastSupportTimes) < fastS + 1
        error("GriddedStreamfunction:InsufficientFastSupportTimes", ...
            "The representative times must contain at least fastS + 1 unique samples.");
    end
    fastKnotPoints = BSpline.knotPointsForDataPoints(fastSupportTimes, S=fastS, splineDOF=numel(fastSupportTimes));
else
    representativeTimes = zeros(0, 1);
    fastKnotPoints = GriddedStreamfunction.validateFastKnotPoints(fastKnotPoints, fastS, allT);
    fastSupportTimes = fitSupportTimes;
end

fastBasisMatrix = BSpline.matrixForDataPoints(allT, knotPoints=fastKnotPoints, S=fastS);
fastBasisMatrix = reshape(fastBasisMatrix, numel(allT), []);
if fastS > 0
    fastDerivativeTensor = BSpline.matrixForDataPoints(allT, knotPoints=fastKnotPoints, S=fastS, D=1);
    fastDerivativeMatrix = reshape(fastDerivativeTensor(:, :, 2), numel(allT), []);
else
    fastDerivativeMatrix = zeros(size(fastBasisMatrix));
end

centerXCoefficients = fastBasisMatrix \ allX;
centerYCoefficients = fastBasisMatrix \ allY;
centerXSpline = TensorSpline(S=fastS, knotPoints=fastKnotPoints, xi=centerXCoefficients);
centerYSpline = TensorSpline(S=fastS, knotPoints=fastKnotPoints, xi=centerYCoefficients);
centerOfMassTrajectory = TrajectorySpline.fromComponentSplines(fitSupportTimes, centerXSpline, centerYSpline);

mxAll = centerOfMassTrajectory.x(allT);
myAll = centerOfMassTrajectory.y(allT);
mxDotAll = fastDerivativeMatrix * centerXCoefficients;
myDotAll = fastDerivativeMatrix * centerYCoefficients;
qAll = allX - mxAll;
rAll = allY - myAll;

if isempty(psiKnotPoints)
    psiKnotPoints = GriddedStreamfunction.defaultPsiKnotPoints(qAll, rAll, allT, psiS);
else
    psiKnotPoints = GriddedStreamfunction.validatePsiKnotPoints(psiKnotPoints);
end

try
    psiBasisSize = TensorSpline.basisSizeFromKnotCell(psiKnotPoints, psiS + 1);
catch
    error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
        "psiKnotPoints must define at least one basis function in each dimension.");
end

qDomain = [psiKnotPoints{1}(1), psiKnotPoints{1}(end)];
rDomain = [psiKnotPoints{2}(1), psiKnotPoints{2}(end)];
tDomain = [psiKnotPoints{3}(1), psiKnotPoints{3}(end)];
if any(qAll < qDomain(1) | qAll > qDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Centered q observations must lie inside the supplied psi spline domain.");
end
if any(rAll < rDomain(1) | rAll > rDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Centered r observations must lie inside the supplied psi spline domain.");
end
if any(allT < tDomain(1) | allT > tDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Observation times must lie inside the supplied psi time domain.");
end

pointMatrix = [qAll, rAll, allT];
uDesignMatrix = -sparse(TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[0 1 0]));
vDesignMatrix = sparse(TensorSpline.matrixForPointMatrix(pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[1 0 0]));
psiDesignMatrix = [uDesignMatrix; vDesignMatrix];
observedCenteredVelocity = [allXDot - mxDotAll; allYDot - myDotAll];

numSpatialCoefficients = prod(psiBasisSize(1:2));
numTimeCoefficients = psiBasisSize(3);
if numSpatialCoefficients == 1
    reducedSpatialBasis = sparse(1, 0);
    reducedBasisMatrix = sparse(prod(psiBasisSize), 0);
    reducedDesignMatrix = sparse(size(psiDesignMatrix, 1), 0);
    reducedCoefficients = zeros(0, 1);
    streamfunctionCoefficients = zeros(prod(psiBasisSize), 1);
else
    reducedSpatialBasis = spalloc(numSpatialCoefficients, numSpatialCoefficients - 1, ...
        2 * (numSpatialCoefficients - 1));
    reducedSpatialBasis(1:(numSpatialCoefficients - 1), :) = speye(numSpatialCoefficients - 1);
    reducedSpatialBasis(numSpatialCoefficients, :) = -1;
    reducedBasisMatrix = kron(speye(numTimeCoefficients), reducedSpatialBasis);
    reducedDesignMatrix = psiDesignMatrix * reducedBasisMatrix;
    reducedCoefficients = reducedDesignMatrix \ observedCenteredVelocity;
    streamfunctionCoefficients = reducedBasisMatrix * reducedCoefficients;
end

streamfunctionSpline = TensorSpline(S=psiS, knotPoints=psiKnotPoints, xi=streamfunctionCoefficients);
uMesoscaleObserved = -streamfunctionSpline.valueAtPoints(qAll, rAll, allT, D=[0 1 0]);
vMesoscaleObserved = streamfunctionSpline.valueAtPoints(qAll, rAll, allT, D=[1 0 0]);
rhoX = allXDot - uMesoscaleObserved;
rhoY = allYDot - vMesoscaleObserved;

backgroundXCoefficients = fastBasisMatrix \ rhoX;
backgroundYCoefficients = fastBasisMatrix \ rhoY;
backgroundXSpline = TensorSpline(S=fastS, knotPoints=fastKnotPoints, xi=backgroundXCoefficients);
backgroundYSpline = TensorSpline(S=fastS, knotPoints=fastKnotPoints, xi=backgroundYCoefficients);
backgroundVelocityTrajectory = TrajectorySpline.fromComponentSplines(fitSupportTimes, backgroundXSpline, backgroundYSpline);

uBackgroundObserved = backgroundVelocityTrajectory.x(allT);
vBackgroundObserved = backgroundVelocityTrajectory.y(allT);
submesoscaleX = rhoX - uBackgroundObserved;
submesoscaleY = rhoY - vBackgroundObserved;

self.streamfunctionSpline = streamfunctionSpline;
self.centerOfMassTrajectory = centerOfMassTrajectory;
self.backgroundVelocityTrajectory = backgroundVelocityTrajectory;
self.fastKnotPoints = fastKnotPoints;
self.fastS = fastS;
self.psiKnotPoints = psiKnotPoints;
self.psiS = psiS;
self.representativeTimes = representativeTimes;
self.fitSupportTimes = fitSupportTimes;
self.sampleCounts = sampleCounts;
self.observationTimes = allT;
self.observedX = allX;
self.observedY = allY;
self.observedXVelocity = allXDot;
self.observedYVelocity = allYDot;
self.centeredX = qAll;
self.centeredY = rAll;
self.centerVelocityX = mxDotAll;
self.centerVelocityY = myDotAll;
self.uMesoscaleObserved = uMesoscaleObserved;
self.vMesoscaleObserved = vMesoscaleObserved;
self.rhoX = rhoX;
self.rhoY = rhoY;
self.uBackgroundObserved = uBackgroundObserved;
self.vBackgroundObserved = vBackgroundObserved;
self.submesoscaleX = submesoscaleX;
self.submesoscaleY = submesoscaleY;

self.fitDiagnostics = struct( ...
    "trajectories", trajectories, ...
    "fastSupportTimes", fastSupportTimes, ...
    "fastBasisMatrix", fastBasisMatrix, ...
    "fastDerivativeMatrix", fastDerivativeMatrix, ...
    "centerXCoefficients", centerXCoefficients, ...
    "centerYCoefficients", centerYCoefficients, ...
    "psiBasisSize", psiBasisSize, ...
    "psiReducedSpatialBasis", reducedSpatialBasis, ...
    "psiReducedBasisMatrix", reducedBasisMatrix, ...
    "psiReducedDesignMatrix", reducedDesignMatrix, ...
    "psiDesignMatrix", psiDesignMatrix, ...
    "psiObservedCenteredVelocity", observedCenteredVelocity, ...
    "psiReducedCoefficients", reducedCoefficients, ...
    "backgroundXCoefficients", backgroundXCoefficients, ...
    "backgroundYCoefficients", backgroundYCoefficients);
end
