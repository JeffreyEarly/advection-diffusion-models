function fastKnotPoints = validateFastKnotPoints(fastKnotPoints, fastS, allT)
validateattributes(fastKnotPoints, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
fastKnotPoints = reshape(fastKnotPoints, [], 1);
if any(diff(fastKnotPoints) < 0)
    error("GriddedStreamfunction:InvalidFastKnotPoints", ...
        "fastKnotPoints must be nondecreasing.");
end
if numel(fastKnotPoints) <= fastS + 1
    error("GriddedStreamfunction:InvalidFastKnotPoints", ...
        "fastKnotPoints must define at least one fast temporal basis function.");
end
if any(allT < fastKnotPoints(1) | allT > fastKnotPoints(end))
    error("GriddedStreamfunction:ObservationOutsideFastDomain", ...
        "Observation times must lie inside the supplied fast spline domain.");
end
end
