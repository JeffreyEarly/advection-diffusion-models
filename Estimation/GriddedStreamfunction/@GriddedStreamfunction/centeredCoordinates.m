function [tEval, q, r] = centeredCoordinates(self, t, x, y)
if ~isequal(size(x), size(y))
    error("GriddedStreamfunction:EvaluationSizeMismatch", ...
        "x and y must have the same size.");
end

if isscalar(t)
    tEval = repmat(t, size(x));
elseif isvector(t) && ismatrix(x) && size(x, 1) == numel(t)
    tEval = repmat(reshape(t, [], 1), 1, size(x, 2));
else
    if ~isequal(size(t), size(x))
        error("GriddedStreamfunction:EvaluationTimeSizeMismatch", ...
            "t must be scalar or have the same size as x and y.");
    end
    tEval = t;
end

q = x - self.centerOfMassTrajectory.x(tEval);
r = y - self.centerOfMassTrajectory.y(tEval);
end
