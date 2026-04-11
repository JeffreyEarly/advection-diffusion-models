function sampleData = decompositionSampleData(self, trajectories)
arguments (Input)
    self (1,1) GriddedStreamfunction
    trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
end
arguments (Output)
    sampleData struct
end

trajectories = reshape(trajectories, [], 1);
[tCell, xCell, yCell, xDotCell, yDotCell] = trajectorySampleData(trajectories);

allT = vertcat(tCell{:});
allX = vertcat(xCell{:});
allY = vertcat(yCell{:});
allXDot = vertcat(xDotCell{:});
allYDot = vertcat(yDotCell{:});
sampleSupportTimes = unique(allT, "sorted");

validateFastTimeDomain(self, allT);

mxDotAll = reshape(self.centerOfMassTrajectory.u(allT), [], 1);
myDotAll = reshape(self.centerOfMassTrajectory.v(allT), [], 1);
[~, qAll, rAll] = centeredCoordinates(self, allT, allX, allY);
validatePsiDomain(self, qAll, rAll, allT);

uMesoscaleObserved = reshape(self.uMesoscale(allT, allX, allY), [], 1);
vMesoscaleObserved = reshape(self.vMesoscale(allT, allX, allY), [], 1);

fastState = GriddedStreamfunction.fastBasisState(allT, self.fastKnotPoints, self.fastS, false);
% Reuse the fast-basis factorization in the same least-squares projector B * (B \ A).
mesoscaleComObserved = fastState.B * (fastState.solver \ [uMesoscaleObserved, vMesoscaleObserved]);
uMesoscaleComObserved = mesoscaleComObserved(:, 1);
vMesoscaleComObserved = mesoscaleComObserved(:, 2);
uMesoscaleRelativeObserved = uMesoscaleObserved - uMesoscaleComObserved;
vMesoscaleRelativeObserved = vMesoscaleObserved - vMesoscaleComObserved;

uBackgroundObserved = mxDotAll - uMesoscaleComObserved;
vBackgroundObserved = myDotAll - vMesoscaleComObserved;
uCenteredObserved = allXDot - mxDotAll;
vCenteredObserved = allYDot - myDotAll;
uCenteredSubmesoscaleObserved = uCenteredObserved - uMesoscaleRelativeObserved;
vCenteredSubmesoscaleObserved = vCenteredObserved - vMesoscaleRelativeObserved;
uSubmesoscaleObserved = allXDot - uBackgroundObserved - uMesoscaleObserved;
vSubmesoscaleObserved = allYDot - vBackgroundObserved - vMesoscaleObserved;

sampleData = struct( ...
    "tCell", {tCell}, ...
    "xCell", {xCell}, ...
    "yCell", {yCell}, ...
    "xDotCell", {xDotCell}, ...
    "yDotCell", {yDotCell}, ...
    "allT", allT, ...
    "allX", allX, ...
    "allY", allY, ...
    "allXDot", allXDot, ...
    "allYDot", allYDot, ...
    "sampleSupportTimes", sampleSupportTimes, ...
    "mxDotAll", mxDotAll, ...
    "myDotAll", myDotAll, ...
    "qAll", qAll, ...
    "rAll", rAll, ...
    "B", fastState.B, ...
    "uMesoscaleObserved", uMesoscaleObserved, ...
    "vMesoscaleObserved", vMesoscaleObserved, ...
    "uMesoscaleComObserved", uMesoscaleComObserved, ...
    "vMesoscaleComObserved", vMesoscaleComObserved, ...
    "uMesoscaleRelativeObserved", uMesoscaleRelativeObserved, ...
    "vMesoscaleRelativeObserved", vMesoscaleRelativeObserved, ...
    "uBackgroundObserved", uBackgroundObserved, ...
    "vBackgroundObserved", vBackgroundObserved, ...
    "uCenteredObserved", uCenteredObserved, ...
    "vCenteredObserved", vCenteredObserved, ...
    "uCenteredSubmesoscaleObserved", uCenteredSubmesoscaleObserved, ...
    "vCenteredSubmesoscaleObserved", vCenteredSubmesoscaleObserved, ...
    "uSubmesoscaleObserved", uSubmesoscaleObserved, ...
    "vSubmesoscaleObserved", vSubmesoscaleObserved);
end

function [tCell, xCell, yCell, xDotCell, yDotCell] = trajectorySampleData(trajectories)
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
end

function validateFastTimeDomain(self, allT)
if any(allT < self.fastKnotPoints(1) | allT > self.fastKnotPoints(end))
    error("GriddedStreamfunction:ObservationOutsideFastDomain", ...
        "Observation times must lie inside the fitted fast spline domain.");
end
end

function validatePsiDomain(self, qAll, rAll, allT)
qDomain = [self.psiKnotPoints{1}(1), self.psiKnotPoints{1}(end)];
rDomain = [self.psiKnotPoints{2}(1), self.psiKnotPoints{2}(end)];
tDomain = [self.psiKnotPoints{3}(1), self.psiKnotPoints{3}(end)];

if any(qAll < qDomain(1) | qAll > qDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Centered q observations must lie inside the fitted psi spline domain.");
end
if any(rAll < rDomain(1) | rAll > rDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Centered r observations must lie inside the fitted psi spline domain.");
end
if any(allT < tDomain(1) | allT > tDomain(2))
    error("GriddedStreamfunction:ObservationOutsidePsiDomain", ...
        "Observation times must lie inside the fitted psi time domain.");
end
end
