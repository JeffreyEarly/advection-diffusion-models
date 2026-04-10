function kappaEstimate = kappaEstimateFromFit(fit)
sampleData = fit.decompositionSampleData(fit.observedTrajectories);
nTrajectories = numel(sampleData.tCell);
kappaByTrajectory = zeros(nTrajectories, 1);

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
    xEnd = terminalDisplacement(ti, sampleData.uSubmesoscaleObserved(sampleIndices), componentS);
    yEnd = terminalDisplacement(ti, sampleData.vSubmesoscaleObserved(sampleIndices), componentS);
    kappaByTrajectory(iTrajectory) = (xEnd.^2 + yEnd.^2) / (4 * duration);
    sampleStartIndex = sampleStartIndex + nSamples;
end

kappaEstimate = mean(kappaByTrajectory);
end

function displacement = terminalDisplacement(t, velocitySamples, componentS)
velocitySpline = InterpolatingSpline.fromGriddedValues(t, velocitySamples, S=componentS);
integralSpline = cumsum(velocitySpline);
displacement = integralSpline(t(end));
end
