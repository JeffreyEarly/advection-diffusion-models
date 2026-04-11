function filename = caseStudyBootstrapCacheFilename(siteNumber, options)
arguments
    siteNumber (1,1) double {mustBeMember(siteNumber, [1 2])}
    options.nBootstraps (1,1) double {mustBeInteger, mustBePositive}
    options.randomSeed (1,1) double {mustBeInteger, mustBeFinite}
    options.scoreStride (1,1) double {mustBeInteger, mustBePositive}
    options.psiS (1,3) double {mustBeInteger, mustBeNonnegative}
    options.fastS (1,1) double {mustBeInteger, mustBeNonnegative}
    options.mesoscaleConstraint {mustBeTextScalar, mustBeMember(options.mesoscaleConstraint, ["none", "zeroVorticity", "zeroStrain"])}
end

psiSTag = join(string(options.psiS), "-");
filename = "Rho" + siteNumber + ...
    "GriddedStreamfunctionBootstrapFits" + options.nBootstraps + ...
    "_seed" + options.randomSeed + ...
    "_stride" + options.scoreStride + ...
    "_fastS" + options.fastS + ...
    "_psiS" + psiSTag + ...
    "_mesoscaleConstraint-" + string(options.mesoscaleConstraint) + ".mat";
end
