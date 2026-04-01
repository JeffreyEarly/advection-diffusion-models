function representativeTimes = representativeObservationTimes(tCell)
pooledTimes = sort(vertcat(tCell{:}));
representativeTimes = zeros(0, 1);
pooledIndex = 1;

while pooledIndex <= numel(pooledTimes)
    tNow = pooledTimes(pooledIndex);
    representativeTimes(end + 1, 1) = tNow; %#ok<AGROW>
    deployedCount = 0;
    for iTrajectory = 1:numel(tCell)
        ti = tCell{iTrajectory};
        deployedCount = deployedCount + double(tNow >= ti(1) && tNow <= ti(end));
    end
    if deployedCount < 1
        error("GriddedStreamfunction:InvalidRepresentativeTimes", ...
            "The representative-time stride rule produced a zero active-drifter count.");
    end
    pooledIndex = pooledIndex + deployedCount;
end
end
