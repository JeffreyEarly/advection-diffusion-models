function kappaEstimate = kappaEstimateFromFit(fit, terminalWeightsByTrajectory)
if nargin < 2
    terminalWeightsByTrajectory = {};
end

sampleData = fit.decompositionSampleData(fit.observedTrajectories);
nTrajectories = numel(sampleData.tCell);
kappaByTrajectory = zeros(nTrajectories, 1);
if isempty(terminalWeightsByTrajectory)
    terminalWeightsByTrajectory = GriddedStreamfunctionBootstrap.terminalDisplacementWeightsForSampleTimes(sampleData.tCell);
elseif numel(terminalWeightsByTrajectory) ~= nTrajectories
    error("GriddedStreamfunctionBootstrap:InvalidTerminalWeightCount", ...
        "Expected one terminal-weight row for each sampled trajectory.");
end

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
    terminalWeights = reshape(terminalWeightsByTrajectory{iTrajectory}, 1, []);
    xEnd = terminalWeights * reshape(sampleData.uSubmesoscaleObserved(sampleIndices), [], 1);
    yEnd = terminalWeights * reshape(sampleData.vSubmesoscaleObserved(sampleIndices), [], 1);
    kappaByTrajectory(iTrajectory) = (xEnd.^2 + yEnd.^2) / (4 * duration);
    sampleStartIndex = sampleStartIndex + nSamples;
end

kappaEstimate = mean(kappaByTrajectory);
end
