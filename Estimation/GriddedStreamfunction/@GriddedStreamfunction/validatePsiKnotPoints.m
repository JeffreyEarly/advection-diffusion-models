function psiKnotPoints = validatePsiKnotPoints(psiKnotPoints)
if ~iscell(psiKnotPoints) || numel(psiKnotPoints) ~= 3
    error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
        "psiKnotPoints must be a cell array {qKnot, rKnot, tKnot}.");
end

psiKnotPoints = reshape(psiKnotPoints, 1, []);
for iDim = 1:3
    validateattributes(psiKnotPoints{iDim}, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
    psiKnotPoints{iDim} = reshape(psiKnotPoints{iDim}, [], 1);
    if any(diff(psiKnotPoints{iDim}) < 0)
        error("GriddedStreamfunction:InvalidPsiKnotPoints", ...
            "Each psi knot vector must be nondecreasing.");
    end
end
end
