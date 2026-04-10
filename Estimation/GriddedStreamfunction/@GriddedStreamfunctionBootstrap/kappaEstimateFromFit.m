function kappaEstimate = kappaEstimateFromFit(fit)
sampleData = fit.decompositionSampleData(fit.observedTrajectories);
nTrajectories = numel(sampleData.tCell);
kappaByTrajectory = zeros(nTrajectories, 1);
cachedTimes = cell(0, 1);
cachedDegrees = zeros(0, 1);
cachedTerminalWeights = cell(0, 1);

sampleStartIndex = 1;
for iTrajectory = 1:nTrajectories
    ti = sampleData.tCell{iTrajectory};
    duration = ti(end) - ti(1);
    if duration <= 0
        error("GriddedStreamfunctionBootstrap:InvalidTrajectoryDuration", ...
            "Each submesoscale trajectory must span a strictly positive time interval.");
    end

    nSamples = numel(ti);
    sampleIndices = sampleStartIndex:(sampleStartIndex + nSamples - 1);
    componentS = min(3, nSamples - 1);
    [terminalWeights, cachedTimes, cachedDegrees, cachedTerminalWeights] = terminalDisplacementWeights( ...
        ti, componentS, cachedTimes, cachedDegrees, cachedTerminalWeights);
    xEnd = terminalWeights * reshape(sampleData.uSubmesoscaleObserved(sampleIndices), [], 1);
    yEnd = terminalWeights * reshape(sampleData.vSubmesoscaleObserved(sampleIndices), [], 1);
    kappaByTrajectory(iTrajectory) = (xEnd.^2 + yEnd.^2) / (4 * duration);
    sampleStartIndex = sampleStartIndex + nSamples;
end

kappaEstimate = mean(kappaByTrajectory);
end

function [terminalWeights, cachedTimes, cachedDegrees, cachedTerminalWeights] = terminalDisplacementWeights( ...
        t, componentS, cachedTimes, cachedDegrees, cachedTerminalWeights)
t = reshape(t, [], 1);
for iCache = 1:numel(cachedTimes)
    if cachedDegrees(iCache) == componentS && isequal(cachedTimes{iCache}, t)
        terminalWeights = cachedTerminalWeights{iCache};
        return
    end
end

terminalWeights = computeTerminalDisplacementWeights(t, componentS);
cachedTimes{end + 1, 1} = t;
cachedDegrees(end + 1, 1) = componentS;
cachedTerminalWeights{end + 1, 1} = terminalWeights;
end

function terminalWeights = computeTerminalDisplacementWeights(t, componentS)
tKnot = BSpline.knotPointsForDataPoints(t, S=componentS);
terminalWeights = BSpline.integralMatrixForDataPoints(t, t(end), knotPoints=tKnot, S=componentS);
end
