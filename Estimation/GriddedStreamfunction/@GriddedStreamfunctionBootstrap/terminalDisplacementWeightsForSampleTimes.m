function terminalWeightsByTrajectory = terminalDisplacementWeightsForSampleTimes(tCell)
nTrajectories = numel(tCell);
terminalWeightsByTrajectory = cell(nTrajectories, 1);
cachedTimes = cell(0, 1);
cachedDegrees = zeros(0, 1);
cachedTerminalWeights = cell(0, 1);

for iTrajectory = 1:nTrajectories
    ti = reshape(tCell{iTrajectory}, [], 1);
    componentS = min(3, numel(ti) - 1);
    [terminalWeights, cachedTimes, cachedDegrees, cachedTerminalWeights] = terminalDisplacementWeights( ...
        ti, componentS, cachedTimes, cachedDegrees, cachedTerminalWeights);
    terminalWeightsByTrajectory{iTrajectory} = terminalWeights;
end
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
