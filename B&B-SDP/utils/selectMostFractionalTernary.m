function [idx, fracval] = selectMostFractionalTernary(x)

%   Choose the index of the most fractional variable
%   [idx, fracval] = selectMostFractionalTernary(x)
%   - x : relaxation solution vector (assumed in [-1,1])
%
%   Returns:
%   - idx     : index of the most fractional variable (empty if all integral)
%   - fracval : fractionality value at that index
%
%   Fractionality = distance to nearest admissible integer in {-1,0,1}.
%   The variable with the largest fractionality is selected.

    x = x(:); % ensure column vector

    % Identify integral ternary values
    isInt = isIntegralTernary(x);

    if all(isInt)
        idx = -1;
        fracval = 0;
        return;
    end

    % Distance to nearest admissible integer in {-1,0,1}
    dminus1 = abs(x + 1);
    dzero   = abs(x);
    done    = abs(x - 1);
    frac    = min([dminus1, dzero, done], [], 2);

    % Ignore integral variables
    frac(isInt) = -inf;

    % Select variable with largest fractionality
    [fracval, idx] = max(frac);

    if ~isfinite(fracval) || fracval <= 0
        idx = [];
        fracval = 0;
    end
end

function t = isIntegralTernary(x)
% True if x is (numerically) in {-1,0,1}
    tol = 1e-8;
    t = abs(x + 1) <= tol | abs(x) <= tol | abs(x - 1) <= tol;
end
