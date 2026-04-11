function bootstrap = loadOrBuildCaseStudyBootstrap(siteNumber, trajectories, options)
arguments
    siteNumber (1,1) double {mustBeMember(siteNumber, [1 2])}
    trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
    options.nBootstraps (1,1) double {mustBeInteger, mustBePositive} = 1000
    options.randomSeed (1,1) double {mustBeInteger, mustBeFinite} = 0
    options.scoreStride (1,1) double {mustBeInteger, mustBePositive} = 6
    options.psiS (1,3) double {mustBeInteger, mustBeNonnegative} = [2 2 4]
    options.fastS (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    options.mesoscaleConstraint {mustBeTextScalar, mustBeMember(options.mesoscaleConstraint, ["none", "zeroVorticity", "zeroStrain"])} = "none"
end

trajectories = reshape(trajectories, [], 1);
scriptDir = fileparts(mfilename("fullpath"));
bootstrapDataDir = fullfile(scriptDir, "BootstrapData");
if exist(bootstrapDataDir, "dir") == 0
    mkdir(bootstrapDataDir);
end

cachePath = fullfile( ...
    bootstrapDataDir, ...
    caseStudyBootstrapCacheFilename( ...
        siteNumber, ...
        nBootstraps=options.nBootstraps, ...
        randomSeed=options.randomSeed, ...
        scoreStride=options.scoreStride, ...
        psiS=options.psiS, ...
        fastS=options.fastS, ...
        mesoscaleConstraint=options.mesoscaleConstraint));

if isfile(cachePath)
    try
        fileData = load(cachePath, "bootstrap");
        if isfield(fileData, "bootstrap") && isValidBootstrap(fileData.bootstrap, trajectories, options)
            bootstrap = fileData.bootstrap;
            return
        end
    catch
    end
end

bootstrap = GriddedStreamfunctionBootstrap( ...
    trajectories, ...
    nBootstraps=options.nBootstraps, ...
    randomSeed=options.randomSeed, ...
    scoreStride=options.scoreStride, ...
    psiS=options.psiS, ...
    fastS=options.fastS, ...
    mesoscaleConstraint=options.mesoscaleConstraint);
save(cachePath, "bootstrap");
end

function tf = isValidBootstrap(bootstrap, trajectories, options)
tf = isa(bootstrap, "GriddedStreamfunctionBootstrap");
if ~tf
    return
end

expectedQueryTimes = defaultQueryTimes(trajectories);
expectedScoreTimes = expectedQueryTimes(1:options.scoreStride:end);

tf = bootstrap.nBootstraps == options.nBootstraps && ...
    bootstrap.randomSeed == options.randomSeed && ...
    isequal(bootstrap.fullFit.psiS, options.psiS) && ...
    bootstrap.fullFit.fastS == options.fastS && ...
    bootstrap.fullFit.mesoscaleConstraint == string(options.mesoscaleConstraint) && ...
    sameVector(bootstrap.queryTimes, expectedQueryTimes) && ...
    sameVector(bootstrap.scoreTimes, expectedScoreTimes) && ...
    sameTrajectories(bootstrap.observedTrajectories, trajectories);
end

function queryTimes = defaultQueryTimes(trajectories)
nTrajectories = numel(trajectories);
allTimes = cell(nTrajectories, 1);

for iTrajectory = 1:nTrajectories
    allTimes{iTrajectory} = reshape(trajectories(iTrajectory).t, [], 1);
end

pooledTimes = unique(vertcat(allTimes{:}), "sorted");
overlapStart = max(cellfun(@(ti) ti(1), allTimes));
overlapEnd = min(cellfun(@(ti) ti(end), allTimes));
queryTimes = pooledTimes(pooledTimes >= overlapStart & pooledTimes <= overlapEnd);
end

function tf = sameTrajectories(referenceTrajectories, trajectories)
if numel(referenceTrajectories) ~= numel(trajectories)
    tf = false;
    return
end

tf = true;
for iTrajectory = 1:numel(trajectories)
    ti = reshape(trajectories(iTrajectory).t, [], 1);
    referenceTimes = reshape(referenceTrajectories(iTrajectory).t, [], 1);
    if ~sameVector(referenceTimes, ti)
        tf = false;
        return
    end

    if ~sameVector(referenceTrajectories(iTrajectory).x(ti), trajectories(iTrajectory).x(ti)) || ...
            ~sameVector(referenceTrajectories(iTrajectory).y(ti), trajectories(iTrajectory).y(ti))
        tf = false;
        return
    end
end
end

function tf = sameVector(a, b)
a = reshape(a, [], 1);
b = reshape(b, [], 1);
tf = isequal(size(a), size(b)) && max(abs(a - b), [], "all") <= 1e-12;
end
