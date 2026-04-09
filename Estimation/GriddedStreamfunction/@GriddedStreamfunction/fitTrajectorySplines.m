function fitTrajectorySplines(self, trajectories, psiKnotPoints, psiS, fastKnotPoints, fastS)
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
fitSupportTimes = unique(allT, "sorted");

if isempty(fastKnotPoints)
    representativeTimes = GriddedStreamfunction.representativeObservationTimes(tCell);
    fastSupportTimes = unique(representativeTimes, "sorted");
    if numel(fastSupportTimes) < fastS + 1
        error("GriddedStreamfunction:InsufficientFastSupportTimes", ...
            "The representative times must contain at least fastS + 1 unique samples.");
    end

    fastKnotPoints = BSpline.knotPointsForDataPoints( ...
        fastSupportTimes, S=fastS, splineDOF=numel(fastSupportTimes));
else
    representativeTimes = zeros(0, 1);
    if any(allT < fastKnotPoints(1) | allT > fastKnotPoints(end))
        error("GriddedStreamfunction:ObservationOutsideFastDomain", ...
            "Observation times must lie inside the supplied fast spline domain.");
    end
end

B = reshape(BSpline.matrixForDataPoints(allT, knotPoints=fastKnotPoints, S=fastS), ...
    numel(allT), []);
if fastS > 0
    dBTensor = BSpline.matrixForDataPoints(allT, knotPoints=fastKnotPoints, S=fastS, D=1);
    dB = reshape(dBTensor(:, :, 2), numel(allT), []);
else
    dB = zeros(size(B));
end

bX = B \ allX;
bY = B \ allY;
centerXSpline = TensorSpline(S=fastS, knotPoints=fastKnotPoints, xi=bX);
centerYSpline = TensorSpline(S=fastS, knotPoints=fastKnotPoints, xi=bY);
centerOfMassTrajectory = TrajectorySpline.fromComponentSplines(fitSupportTimes, centerXSpline, centerYSpline);

mxAll = B * bX;
myAll = B * bY;
% \dot{m} comes from differentiating the fitted fast COM spline.
mxDotAll = dB * bX;
myDotAll = dB * bY;
qAll = allX - mxAll;
rAll = allY - myAll;

if isempty(psiKnotPoints)
    psiKnotPoints = GriddedStreamfunction.defaultPsiKnotPoints(qAll, rAll, allT, psiS);
end

psiBasisSize = TensorSpline.basisSizeFromKnotCell(psiKnotPoints, psiS + 1);
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
Hx = -sparse(TensorSpline.matrixForPointMatrix( ...
    pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[0 1 0]));
Hy = sparse(TensorSpline.matrixForPointMatrix( ...
    pointMatrix, knotPoints=psiKnotPoints, S=psiS, D=[1 0 0]));

numSpatialCoefficients = prod(psiBasisSize(1:2));
if numSpatialCoefficients == 1
    Gspatial = sparse(1, 0);
else
    spatialSupportCounts = max(2 * psiBasisSize(1:2), psiBasisSize(1:2) + 2);
    qGaugePoints = linspace(qDomain(1), qDomain(2), spatialSupportCounts(1)).';
    rGaugePoints = linspace(rDomain(1), rDomain(2), spatialSupportCounts(2)).';
    [qGaugeGrid, rGaugeGrid] = ndgrid(qGaugePoints, rGaugePoints);
    spatialPointMatrix = [qGaugeGrid(:), rGaugeGrid(:)];
    spatialBasisMatrix = TensorSpline.matrixForPointMatrix( ...
        spatialPointMatrix, knotPoints=psiKnotPoints(1:2), S=psiS(1:2));
    constantGaugeVector = spatialBasisMatrix \ ones(size(qGaugeGrid(:)));
    constantGaugeResidual = spatialBasisMatrix * constantGaugeVector - ones(size(qGaugeGrid(:)));
    if norm(constantGaugeResidual, inf) > 1e-10
        error("GriddedStreamfunction:InvalidPsiGauge", ...
            "The spatial spline basis must represent a constant streamfunction gauge.");
    end

    Gspatial = sparse(null(reshape(constantGaugeVector, [], 1).', "r"));
end
G = kron(speye(psiBasisSize(3)), Gspatial);

Hx = Hx * G;
Hy = Hy * G;
projectFast = @(A) B * (B \ A); % S(A) = B * (B \ A)
McomX = sparse(projectFast(Hx));
McomY = sparse(projectFast(Hy));
Mcom = [McomX; McomY];
Mrel = [Hx - McomX; Hy - McomY];
dcom = [mxDotAll; myDotAll];
drel = [allXDot - mxDotAll; allYDot - myDotAll];
% d = [Mcom; Mrel] * alpha.
alpha = [Mcom; Mrel] \ [dcom; drel];
streamfunctionCoefficients = G * alpha;

streamfunctionSpline = TensorSpline(S=psiS, knotPoints=psiKnotPoints, ...
    xi=streamfunctionCoefficients);
uMesoscaleObserved = full(Hx * alpha);
vMesoscaleObserved = full(Hy * alpha);
uMesoscaleComObserved = projectFast(uMesoscaleObserved);
vMesoscaleComObserved = projectFast(vMesoscaleObserved);
uMesoscaleRelativeObserved = uMesoscaleObserved - uMesoscaleComObserved;
vMesoscaleRelativeObserved = vMesoscaleObserved - vMesoscaleComObserved;

uBackgroundObserved = mxDotAll - uMesoscaleComObserved;
vBackgroundObserved = myDotAll - vMesoscaleComObserved;
backgroundXCoefficients = B \ uBackgroundObserved;
backgroundYCoefficients = B \ vBackgroundObserved;
backgroundXSpline = TensorSpline(S=fastS, knotPoints=fastKnotPoints, xi=backgroundXCoefficients);
backgroundYSpline = TensorSpline(S=fastS, knotPoints=fastKnotPoints, xi=backgroundYCoefficients);

uCenteredSubmesoscaleObserved = drel(1:numel(allT)) - uMesoscaleRelativeObserved;
vCenteredSubmesoscaleObserved = drel((numel(allT) + 1):end) - vMesoscaleRelativeObserved;
uSubmesoscaleObserved = allXDot - uBackgroundObserved - uMesoscaleObserved;
vSubmesoscaleObserved = allYDot - vBackgroundObserved - vMesoscaleObserved;

backgroundTrajectoryXSpline = cumsum(backgroundXSpline);
backgroundTrajectoryYSpline = cumsum(backgroundYSpline);
backgroundTrajectoryXSpline = TensorSpline(S=backgroundTrajectoryXSpline.S, ...
    knotPoints=backgroundTrajectoryXSpline.knotPoints, xi=backgroundTrajectoryXSpline.xi, ...
    xMean=-backgroundTrajectoryXSpline(fitSupportTimes(1)));
backgroundTrajectoryYSpline = TensorSpline(S=backgroundTrajectoryYSpline.S, ...
    knotPoints=backgroundTrajectoryYSpline.knotPoints, xi=backgroundTrajectoryYSpline.xi, ...
    xMean=-backgroundTrajectoryYSpline(fitSupportTimes(1)));
backgroundTrajectory = TrajectorySpline.fromComponentSplines( ...
    fitSupportTimes, backgroundTrajectoryXSpline, backgroundTrajectoryYSpline);

backgroundTrajectories = trajectories;
mesoscaleTrajectories = trajectories;
submesoscaleTrajectories = trajectories;
centeredMesoscaleTrajectories = trajectories;
centeredSubmesoscaleTrajectories = trajectories;

sampleStartIndex = 1;
for iTrajectory = 1:nTrajectories
    ti = tCell{iTrajectory};
    nSamples = numel(ti);
    sampleIndices = sampleStartIndex:(sampleStartIndex + nSamples - 1);
    componentS = min(3, nSamples - 1);
    x0 = allX(sampleIndices(1));
    y0 = allY(sampleIndices(1));
    q0 = qAll(sampleIndices(1));
    r0 = rAll(sampleIndices(1));

    backgroundTrajectories(iTrajectory) = reanchoredBackgroundTrajectory(backgroundTrajectory, ti);
    mesoscaleTrajectories(iTrajectory) = sampledVelocityTrajectory( ...
        ti, uMesoscaleObserved(sampleIndices), vMesoscaleObserved(sampleIndices), componentS, x0, y0);
    submesoscaleTrajectories(iTrajectory) = sampledVelocityTrajectory( ...
        ti, uSubmesoscaleObserved(sampleIndices), vSubmesoscaleObserved(sampleIndices), componentS, 0, 0);
    centeredMesoscaleTrajectories(iTrajectory) = sampledVelocityTrajectory( ...
        ti, uMesoscaleRelativeObserved(sampleIndices), vMesoscaleRelativeObserved(sampleIndices), componentS, q0, r0);
    centeredSubmesoscaleTrajectories(iTrajectory) = sampledVelocityTrajectory( ...
        ti, uCenteredSubmesoscaleObserved(sampleIndices), vCenteredSubmesoscaleObserved(sampleIndices), componentS, 0, 0);

    sampleStartIndex = sampleStartIndex + nSamples;
end

selfDecomposition = struct( ...
    "fixedFrame", struct( ...
        "background", backgroundTrajectories, ...
        "mesoscale", mesoscaleTrajectories, ...
        "submesoscale", submesoscaleTrajectories), ...
    "centeredFrame", struct( ...
        "mesoscale", centeredMesoscaleTrajectories, ...
        "submesoscale", centeredSubmesoscaleTrajectories));

self.streamfunctionSpline = streamfunctionSpline;
self.observedTrajectories = trajectories;
self.centerOfMassTrajectory = centerOfMassTrajectory;
self.backgroundTrajectory = backgroundTrajectory;
self.decomposition = selfDecomposition;
self.fastKnotPoints = fastKnotPoints;
self.fastS = fastS;
self.psiKnotPoints = psiKnotPoints;
self.psiS = psiS;
self.representativeTimes = representativeTimes;
self.fitSupportTimes = fitSupportTimes;
end

function trajectory = sampledVelocityTrajectory(t, uSamples, vSamples, componentS, xMean, yMean)
xVelocitySpline = InterpolatingSpline(t, uSamples, S=componentS);
yVelocitySpline = InterpolatingSpline(t, vSamples, S=componentS);
xIntegralSpline = cumsum(xVelocitySpline);
yIntegralSpline = cumsum(yVelocitySpline);
xPathSpline = TensorSpline(S=xIntegralSpline.S, knotPoints=xIntegralSpline.knotPoints, ...
    xi=xIntegralSpline.xi, xMean=xMean);
yPathSpline = TensorSpline(S=yIntegralSpline.S, knotPoints=yIntegralSpline.knotPoints, ...
    xi=yIntegralSpline.xi, xMean=yMean);
trajectory = TrajectorySpline.fromComponentSplines(t, xPathSpline, yPathSpline);
end

function trajectory = reanchoredBackgroundTrajectory(backgroundTrajectory, ti)
xSpline = TensorSpline(S=backgroundTrajectory.x.S, knotPoints=backgroundTrajectory.x.knotPoints, ...
    xi=backgroundTrajectory.x.xi, xMean=-backgroundTrajectory.x(ti(1)));
ySpline = TensorSpline(S=backgroundTrajectory.y.S, knotPoints=backgroundTrajectory.y.knotPoints, ...
    xi=backgroundTrajectory.y.xi, xMean=-backgroundTrajectory.y(ti(1)));
trajectory = TrajectorySpline.fromComponentSplines(ti, xSpline, ySpline);
end
