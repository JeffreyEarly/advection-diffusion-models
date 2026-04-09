function kappaEstimate = kappaEstimateFromFit(fit)
nTrajectories = numel(fit.decomposition.fixedFrame.submesoscale);
kappaByTrajectory = zeros(nTrajectories, 1);

for iTrajectory = 1:nTrajectories
    trajectory = fit.decomposition.fixedFrame.submesoscale(iTrajectory);
    ti = reshape(trajectory.t, [], 1);
    duration = ti(end) - ti(1);
    if duration <= 0
        error("GriddedStreamfunctionBootstrap:InvalidTrajectoryDuration", ...
            "Each submesoscale trajectory must span a strictly positive time interval.");
    end

    xEnd = trajectory.x(ti(end));
    yEnd = trajectory.y(ti(end));
    kappaByTrajectory(iTrajectory) = (xEnd.^2 + yEnd.^2) / (4 * duration);
end

kappaEstimate = mean(kappaByTrajectory);
end
