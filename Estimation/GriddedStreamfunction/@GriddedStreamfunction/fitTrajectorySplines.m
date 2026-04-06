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
    validateattributes(fastKnotPoints, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
    fastKnotPoints = reshape(fastKnotPoints, [], 1);
    if any(diff(fastKnotPoints) < 0)
        error("GriddedStreamfunction:InvalidFastKnotPoints", ...
            "fastKnotPoints must be nondecreasing.");
    end
    if numel(fastKnotPoints) <= fastS + 1
        error("GriddedStreamfunction:InvalidFastKnotPoints", ...
            "fastKnotPoints must define at least one fast temporal basis function.");
    end
    if any(allT < fastKnotPoints(1) | allT > fastKnotPoints(end))
        error("GriddedStreamfunction:ObservationOutsideFastDomain", ...
            "Observation times must lie inside the supplied fast spline domain.");
    end

    representativeTimes = zeros(0, 1);
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

mxAll = fastBasisMatrix * centerXCoefficients;
myAll = fastBasisMatrix * centerYCoefficients;
mxDotAll = fastDerivativeMatrix * centerXCoefficients;
myDotAll = fastDerivativeMatrix * centerYCoefficients;
qAll = allX - mxAll;
rAll = allY - myAll;

if isempty(psiKnotPoints)
    psiKnotPoints = GriddedStreamfunction.defaultPsiKnotPoints(qAll, rAll, allT, psiS);
else
    if ~iscell(psiKnotPoints) || numel(psiKnotPoints) ~= 3
        error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
            "psiKnotPoints must be a cell array {qKnot, rKnot, tKnot}.");
    end

    psiKnotPoints = reshape(psiKnotPoints, 1, []);
    for iDim = 1:3
        validateattributes(psiKnotPoints{iDim}, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
        psiKnotPoints{iDim} = reshape(psiKnotPoints{iDim}, [], 1);
        if any(diff(psiKnotPoints{iDim}) < 0)
            error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
                "Each psi knot vector must be nondecreasing.");
        end
    end
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

numSpatialCoefficients = prod(psiBasisSize(1:2));
numTimeCoefficients = psiBasisSize(3);
if numSpatialCoefficients == 1
    psiGaugeModeMatrix = 1;
    reducedSpatialBasis = sparse(1, 0);
else
    spatialSupportCounts = max(2 * psiBasisSize(1:2), psiBasisSize(1:2) + 2);
    qGaugePoints = linspace(qDomain(1), qDomain(2), spatialSupportCounts(1)).';
    rGaugePoints = linspace(rDomain(1), rDomain(2), spatialSupportCounts(2)).';
    [qGaugeGrid, rGaugeGrid] = ndgrid(qGaugePoints, rGaugePoints);
    spatialPointMatrix = [qGaugeGrid(:), rGaugeGrid(:)];
    spatialBasisMatrix = TensorSpline.matrixForPointMatrix(spatialPointMatrix, knotPoints=psiKnotPoints(1:2), S=psiS(1:2));
    constantGaugeVector = spatialBasisMatrix \ ones(size(qGaugeGrid(:)));
    constantGaugeResidual = spatialBasisMatrix * constantGaugeVector - ones(size(qGaugeGrid(:)));
    if norm(constantGaugeResidual, inf) > 1e-10
        error("GriddedStreamfunction:InvalidPsiGauge", ...
            "The spatial spline basis must represent a constant streamfunction gauge.");
    end
    psiGaugeModeMatrix = reshape(constantGaugeVector, [], 1);
    reducedSpatialBasis = sparse(null(psiGaugeModeMatrix.', "r"));
end
reducedBasisMatrix = kron(speye(numTimeCoefficients), reducedSpatialBasis);

reducedUDesignMatrix = uDesignMatrix * reducedBasisMatrix;
reducedVDesignMatrix = vDesignMatrix * reducedBasisMatrix;
mesoscaleComXDesignMatrix = sparse(fastBasisMatrix * (fastBasisMatrix \ reducedUDesignMatrix));
mesoscaleComYDesignMatrix = sparse(fastBasisMatrix * (fastBasisMatrix \ reducedVDesignMatrix));
mesoscaleRelativeXDesignMatrix = reducedUDesignMatrix - mesoscaleComXDesignMatrix;
mesoscaleRelativeYDesignMatrix = reducedVDesignMatrix - mesoscaleComYDesignMatrix;

mesoscaleComDesignMatrix = [mesoscaleComXDesignMatrix; mesoscaleComYDesignMatrix];
mesoscaleRelativeDesignMatrix = [mesoscaleRelativeXDesignMatrix; mesoscaleRelativeYDesignMatrix];
comVelocityObserved = [mxDotAll; myDotAll];
relativeXObserved = allXDot - mxDotAll;
relativeYObserved = allYDot - myDotAll;
relativeVelocityObserved = [relativeXObserved; relativeYObserved];
mesoscaleRightHandSide = [comVelocityObserved; relativeVelocityObserved];
mesoscaleDesignMatrix = [mesoscaleComDesignMatrix; mesoscaleRelativeDesignMatrix];
reducedCoefficients = mesoscaleDesignMatrix \ mesoscaleRightHandSide;
streamfunctionCoefficients = reducedBasisMatrix * reducedCoefficients;

streamfunctionSpline = TensorSpline(S=psiS, knotPoints=psiKnotPoints, xi=streamfunctionCoefficients);
uMesoscaleObserved = uDesignMatrix * streamfunctionCoefficients;
vMesoscaleObserved = vDesignMatrix * streamfunctionCoefficients;
uMesoscaleComObserved = fastBasisMatrix * (fastBasisMatrix \ uMesoscaleObserved);
vMesoscaleComObserved = fastBasisMatrix * (fastBasisMatrix \ vMesoscaleObserved);
uMesoscaleRelativeObserved = uMesoscaleObserved - uMesoscaleComObserved;
vMesoscaleRelativeObserved = vMesoscaleObserved - vMesoscaleComObserved;

uBackgroundObserved = mxDotAll - uMesoscaleComObserved;
vBackgroundObserved = myDotAll - vMesoscaleComObserved;
backgroundXCoefficients = fastBasisMatrix \ uBackgroundObserved;
backgroundYCoefficients = fastBasisMatrix \ vBackgroundObserved;
backgroundXSpline = TensorSpline(S=fastS, knotPoints=fastKnotPoints, xi=backgroundXCoefficients);
backgroundYSpline = TensorSpline(S=fastS, knotPoints=fastKnotPoints, xi=backgroundYCoefficients);
backgroundVelocityTrajectory = TrajectorySpline.fromComponentSplines(fitSupportTimes, backgroundXSpline, backgroundYSpline);

uCenteredSubmesoscaleObserved = relativeXObserved - uMesoscaleRelativeObserved;
vCenteredSubmesoscaleObserved = relativeYObserved - vMesoscaleRelativeObserved;
uSubmesoscaleObserved = allXDot - uBackgroundObserved - uMesoscaleObserved;
vSubmesoscaleObserved = allYDot - vBackgroundObserved - vMesoscaleObserved;

backgroundTrajectoryXSpline = cumsum(backgroundXSpline);
backgroundTrajectoryYSpline = cumsum(backgroundYSpline);
backgroundTrajectories = TrajectorySpline.empty(0, 1);
mesoscaleTrajectories = TrajectorySpline.empty(0, 1);
submesoscaleTrajectories = TrajectorySpline.empty(0, 1);
centeredMesoscaleTrajectories = TrajectorySpline.empty(0, 1);
centeredSubmesoscaleTrajectories = TrajectorySpline.empty(0, 1);

sampleStartIndex = 1;
for iTrajectory = 1:nTrajectories
    ti = tCell{iTrajectory};
    nSamples = numel(ti);
    sampleIndices = sampleStartIndex:(sampleStartIndex + nSamples - 1);
    componentS = min(3, numel(ti) - 1);

    backgroundXPathSpline = TensorSpline(S=backgroundTrajectoryXSpline.S, ...
        knotPoints=backgroundTrajectoryXSpline.knotPoints, xi=backgroundTrajectoryXSpline.xi, ...
        xMean=-backgroundTrajectoryXSpline(ti(1)));
    backgroundYPathSpline = TensorSpline(S=backgroundTrajectoryYSpline.S, ...
        knotPoints=backgroundTrajectoryYSpline.knotPoints, xi=backgroundTrajectoryYSpline.xi, ...
        xMean=-backgroundTrajectoryYSpline(ti(1)));
    backgroundTrajectories(end + 1, 1) = TrajectorySpline.fromComponentSplines(ti, backgroundXPathSpline, backgroundYPathSpline); %#ok<AGROW>

    mesoscaleVelocityXSpline = InterpolatingSpline(ti, uMesoscaleObserved(sampleIndices), S=componentS);
    mesoscaleVelocityYSpline = InterpolatingSpline(ti, vMesoscaleObserved(sampleIndices), S=componentS);
    mesoscaleXPathIntegral = cumsum(mesoscaleVelocityXSpline);
    mesoscaleYPathIntegral = cumsum(mesoscaleVelocityYSpline);
    mesoscaleXPathSpline = TensorSpline(S=mesoscaleXPathIntegral.S, ...
        knotPoints=mesoscaleXPathIntegral.knotPoints, xi=mesoscaleXPathIntegral.xi, ...
        xMean=allX(sampleIndices(1)));
    mesoscaleYPathSpline = TensorSpline(S=mesoscaleYPathIntegral.S, ...
        knotPoints=mesoscaleYPathIntegral.knotPoints, xi=mesoscaleYPathIntegral.xi, ...
        xMean=allY(sampleIndices(1)));
    mesoscaleTrajectories(end + 1, 1) = TrajectorySpline.fromComponentSplines(ti, mesoscaleXPathSpline, mesoscaleYPathSpline); %#ok<AGROW>

    submesoscaleVelocityXSpline = InterpolatingSpline(ti, uSubmesoscaleObserved(sampleIndices), S=componentS);
    submesoscaleVelocityYSpline = InterpolatingSpline(ti, vSubmesoscaleObserved(sampleIndices), S=componentS);
    submesoscaleXPathIntegral = cumsum(submesoscaleVelocityXSpline);
    submesoscaleYPathIntegral = cumsum(submesoscaleVelocityYSpline);
    submesoscaleXPathSpline = TensorSpline(S=submesoscaleXPathIntegral.S, ...
        knotPoints=submesoscaleXPathIntegral.knotPoints, xi=submesoscaleXPathIntegral.xi);
    submesoscaleYPathSpline = TensorSpline(S=submesoscaleYPathIntegral.S, ...
        knotPoints=submesoscaleYPathIntegral.knotPoints, xi=submesoscaleYPathIntegral.xi);
    submesoscaleTrajectories(end + 1, 1) = TrajectorySpline.fromComponentSplines(ti, submesoscaleXPathSpline, submesoscaleYPathSpline); %#ok<AGROW>

    centeredMesoscaleVelocityXSpline = InterpolatingSpline(ti, uMesoscaleRelativeObserved(sampleIndices), S=componentS);
    centeredMesoscaleVelocityYSpline = InterpolatingSpline(ti, vMesoscaleRelativeObserved(sampleIndices), S=componentS);
    centeredMesoscaleXPathIntegral = cumsum(centeredMesoscaleVelocityXSpline);
    centeredMesoscaleYPathIntegral = cumsum(centeredMesoscaleVelocityYSpline);
    centeredMesoscaleXPathSpline = TensorSpline(S=centeredMesoscaleXPathIntegral.S, ...
        knotPoints=centeredMesoscaleXPathIntegral.knotPoints, xi=centeredMesoscaleXPathIntegral.xi, ...
        xMean=qAll(sampleIndices(1)));
    centeredMesoscaleYPathSpline = TensorSpline(S=centeredMesoscaleYPathIntegral.S, ...
        knotPoints=centeredMesoscaleYPathIntegral.knotPoints, xi=centeredMesoscaleYPathIntegral.xi, ...
        xMean=rAll(sampleIndices(1)));
    centeredMesoscaleTrajectories(end + 1, 1) = TrajectorySpline.fromComponentSplines(ti, centeredMesoscaleXPathSpline, centeredMesoscaleYPathSpline); %#ok<AGROW>

    centeredSubmesoscaleVelocityXSpline = InterpolatingSpline(ti, uCenteredSubmesoscaleObserved(sampleIndices), S=componentS);
    centeredSubmesoscaleVelocityYSpline = InterpolatingSpline(ti, vCenteredSubmesoscaleObserved(sampleIndices), S=componentS);
    centeredSubmesoscaleXPathIntegral = cumsum(centeredSubmesoscaleVelocityXSpline);
    centeredSubmesoscaleYPathIntegral = cumsum(centeredSubmesoscaleVelocityYSpline);
    centeredSubmesoscaleXPathSpline = TensorSpline(S=centeredSubmesoscaleXPathIntegral.S, ...
        knotPoints=centeredSubmesoscaleXPathIntegral.knotPoints, xi=centeredSubmesoscaleXPathIntegral.xi);
    centeredSubmesoscaleYPathSpline = TensorSpline(S=centeredSubmesoscaleYPathIntegral.S, ...
        knotPoints=centeredSubmesoscaleYPathIntegral.knotPoints, xi=centeredSubmesoscaleYPathIntegral.xi);
    centeredSubmesoscaleTrajectories(end + 1, 1) = TrajectorySpline.fromComponentSplines(ti, centeredSubmesoscaleXPathSpline, centeredSubmesoscaleYPathSpline); %#ok<AGROW>

    sampleStartIndex = sampleStartIndex + nSamples;
end

self.streamfunctionSpline = streamfunctionSpline;
self.centerOfMassTrajectory = centerOfMassTrajectory;
self.backgroundVelocityTrajectory = backgroundVelocityTrajectory;
self.backgroundTrajectories = backgroundTrajectories;
self.mesoscaleTrajectories = mesoscaleTrajectories;
self.submesoscaleTrajectories = submesoscaleTrajectories;
self.centeredMesoscaleTrajectories = centeredMesoscaleTrajectories;
self.centeredSubmesoscaleTrajectories = centeredSubmesoscaleTrajectories;
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
self.uMesoscaleComObserved = uMesoscaleComObserved;
self.vMesoscaleComObserved = vMesoscaleComObserved;
self.uMesoscaleRelativeObserved = uMesoscaleRelativeObserved;
self.vMesoscaleRelativeObserved = vMesoscaleRelativeObserved;
self.uBackgroundObserved = uBackgroundObserved;
self.vBackgroundObserved = vBackgroundObserved;
self.uSubmesoscaleObserved = uSubmesoscaleObserved;
self.vSubmesoscaleObserved = vSubmesoscaleObserved;

self.fitDiagnostics = struct( ...
    "trajectories", trajectories, ...
    "fastSupportTimes", fastSupportTimes, ...
    "fastBasisMatrix", fastBasisMatrix, ...
    "fastDerivativeMatrix", fastDerivativeMatrix, ...
    "centerXCoefficients", centerXCoefficients, ...
    "centerYCoefficients", centerYCoefficients, ...
    "psiBasisSize", psiBasisSize, ...
    "psiGaugeModeMatrix", psiGaugeModeMatrix, ...
    "psiReducedSpatialBasis", reducedSpatialBasis, ...
    "psiReducedBasisMatrix", reducedBasisMatrix, ...
    "uDesignMatrix", uDesignMatrix, ...
    "vDesignMatrix", vDesignMatrix, ...
    "reducedUDesignMatrix", reducedUDesignMatrix, ...
    "reducedVDesignMatrix", reducedVDesignMatrix, ...
    "mesoscaleComXDesignMatrix", mesoscaleComXDesignMatrix, ...
    "mesoscaleComYDesignMatrix", mesoscaleComYDesignMatrix, ...
    "mesoscaleRelativeXDesignMatrix", mesoscaleRelativeXDesignMatrix, ...
    "mesoscaleRelativeYDesignMatrix", mesoscaleRelativeYDesignMatrix, ...
    "mesoscaleComDesignMatrix", mesoscaleComDesignMatrix, ...
    "mesoscaleRelativeDesignMatrix", mesoscaleRelativeDesignMatrix, ...
    "mesoscaleDesignMatrix", mesoscaleDesignMatrix, ...
    "mesoscaleRightHandSide", mesoscaleRightHandSide, ...
    "comVelocityObserved", comVelocityObserved, ...
    "relativeVelocityObserved", relativeVelocityObserved, ...
    "backgroundXCoefficients", backgroundXCoefficients, ...
    "backgroundYCoefficients", backgroundYCoefficients, ...
    "psiReducedCoefficients", reducedCoefficients, ...
    "psiCoefficients", streamfunctionCoefficients, ...
    "uCenteredSubmesoscaleObserved", uCenteredSubmesoscaleObserved, ...
    "vCenteredSubmesoscaleObserved", vCenteredSubmesoscaleObserved);
end
