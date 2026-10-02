function [best_ratio, x_best] = vns_ternary_ratio_qp(A, a, a0, B, b, b0, x0, opts)
% VNS_TERNARY_RATIO_QP
%   Variable Neighborhood Search for
%       min  (x' * A * x + a' * x + a0) / (x' * B * x + b' * x + b0)
%       s.t. x_i in {-1, 0, 1}
%
% INPUT
%   A,B  : (n x n) symmetric matrices
%   a,b  : (n x 1) vectors
%   a0,b0: scalars (constant terms)
%   x0   : (n x 1) initial vector (projected to {-1,0,1})
%   opts : struct with optional fields:
%            .iter_max        (default 3)
%            .kmin            (default 2)
%            .kmax            (default n)
%            .kstep           (default 2)
%            .seed            (default 1)
%            .max_local_iters (default 1e2)
%            .verbose         (default false)
%            .denom_eps       (default 1e-10) % denom must exceed this
%
% OUTPUT
%   x_best     : best solution found (n x 1, entries in {-1,0,1})
%   best_ratio : objective value at x_best

    if nargin < 8
        opts = struct();
    end

    % ----- Basic checks -----
    [n1, n2] = size(A);
    if n1 ~= n2, error('A must be square.'); end
    if any(size(B) ~= [n1,n2]), error('B must be same size as A.'); end
    n = n1;

    if numel(a) ~= n || numel(b) ~= n
        error('Dimension mismatch: a and b must have length size(A,1).');
    end
    if numel(x0) ~= n
        error('Dimension mismatch: x0 must have length size(A,1).');
    end

    A = double(A);  B = double(B);
    a = double(a(:)); b = double(b(:));

    a0 = double(a0);
    b0 = double(b0);
    if ~isscalar(a0) || ~isscalar(b0)
        error('a0 and b0 must be scalars.');
    end

    % ----- Options with defaults -----
    iter_max        = get_opt(opts, 'iter_max',        3);
    kmin            = get_opt(opts, 'kmin',            2);
    kmax            = get_opt(opts, 'kmax',            n);
    kstep           = get_opt(opts, 'kstep',           2);
    seed            = get_opt(opts, 'seed',            1);
    max_local_iters = get_opt(opts, 'max_local_iters', 1e3);
    verbose         = get_opt(opts, 'verbose',         false);
    denom_eps       = get_opt(opts, 'denom_eps',       1e-10);

    % ----- Initial solution and projection to {-1,0,1} -----
    x = double(x0(:));
    x = round(x);
    x(x >  1) =  1;
    x(x < -1) = -1;

    % ----- Initialize aggregates (numerator & denominator) -----
    alphaA = A * x;
    num1   = x' * A * x;   % x'Ax
    num2   = a' * x;       % a'x

    alphaB = B * x;
    den1   = x' * B * x;   % x'Bx
    den2   = b' * x;       % b'x

    best_ratio = safe_ratio(num1 + num2 + a0, den1 + den2 + b0, denom_eps);
    x_best     = x;

    % store best state to restore on non-improvement
    alphaA_b = alphaA; num1_b = num1; num2_b = num2;
    alphaB_b = alphaB; den1_b = den1; den2_b = den2;

    rng(seed);

    % ===================== MAIN VNS LOOP =====================
    for it = 1:iter_max
        k = kmin;
        while k <= kmax

            % ---- Shake ----
            [x, alphaA, num1, num2, alphaB, den1, den2] = shake_ratio( ...
                x, alphaA, num1, num2, alphaB, den1, den2, A, a, B, b, k);

            % ---- Local search in 1-change neighborhood (direct ratio) ----
            [x, alphaA, num1, num2, alphaB, den1, den2] = local_search_ratio( ...
                x, alphaA, num1, num2, alphaB, den1, den2, A, a, a0, B, b, b0, ...
                max_local_iters, denom_eps);

            cur_ratio = safe_ratio(num1 + num2 + a0, den1 + den2 + b0, denom_eps);

            if cur_ratio < best_ratio
                x_best     = x;
                best_ratio = cur_ratio;

                alphaA_b = alphaA; num1_b = num1; num2_b = num2;
                alphaB_b = alphaB; den1_b = den1; den2_b = den2;

                if verbose
                    fprintf('Improved: it=%d k=%d ratio=%.6g\n', it, k, best_ratio);
                end

                k = kmin; % intensify
            else
                % restore best and diversify
                x      = x_best;

                alphaA = alphaA_b; num1 = num1_b; num2 = num2_b;
                alphaB = alphaB_b; den1 = den1_b; den2 = den2_b;

                k = k + kstep;
            end
        end
    end
end

% ======================================================================
% SUBFUNCTIONS
% ======================================================================

function val = get_opt(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = default;
    end
end

function r = safe_ratio(num, den, epsd)
    if den > epsd
        r = num / den;
    else
        r = Inf; % treat invalid/too-small denominator as infeasible/bad
    end
end

% Trial change x_i := t (t in {-1,0,1}, t ~= x_i)
% Updates BOTH numerator and denominator pieces (but does not modify state).
function [num1p, num2p, den1p, den2p] = trial_change_ratio( ...
    i, t, x, num1, num2, den1, den2, alphaA, A, a, alphaB, B, b)

    s = x(i);
    if t == s
        num1p = num1; num2p = num2;
        den1p = den1; den2p = den2;
        return;
    end

    delta = t - s;

    num1p = num1 + 2*delta*alphaA(i) + (delta^2)*A(i,i);
    num2p = num2 + delta*a(i);

    den1p = den1 + 2*delta*alphaB(i) + (delta^2)*B(i,i);
    den2p = den2 + delta*b(i);
end

% Apply change x_i := t (updates x and all aggregates in-place)
function [x, alphaA, num1, num2, alphaB, den1, den2] = apply_change_ratio( ...
    i, t, x, alphaA, num1, num2, A, a, alphaB, den1, den2, B, b)

    s = x(i);
    if t == s
        return;
    end

    delta = t - s;

    % update numerator objective parts
    num1 = num1 + 2*delta*alphaA(i) + (delta^2)*A(i,i);
    num2 = num2 + delta*a(i);

    % update denominator objective parts
    den1 = den1 + 2*delta*alphaB(i) + (delta^2)*B(i,i);
    den2 = den2 + delta*b(i);

    % update alpha vectors
    alphaA = alphaA + delta * A(:,i);
    alphaB = alphaB + delta * B(:,i);

    % update x
    x(i) = t;
end

% Best-improvement local search in 1-change neighborhood (direct ratio)
function [x, alphaA, num1, num2, alphaB, den1, den2] = local_search_ratio( ...
    x, alphaA, num1, num2, alphaB, den1, den2, A, a, a0, B, b, b0, max_iters, denom_eps)

    n   = numel(x);
    EPS = 1e-6;
    it_count = 0;

    improved = true;
    while improved && it_count < max_iters
        it_count = it_count + 1;

        improved   = false;
        best_ratio = safe_ratio(num1 + num2 + a0, den1 + den2 + b0, denom_eps);
        best_i     = 0;
        best_t     = 0;

        for i = 1:n
            s = x(i);
            candidates = [-1,0,1];
            candidates(candidates == s) = [];

            for c_idx = 1:2
                t = candidates(c_idx);

                [num1p, num2p, den1p, den2p] = trial_change_ratio( ...
                    i, t, x, num1, num2, den1, den2, alphaA, A, a, alphaB, B, b);

                r = safe_ratio(num1p + num2p + a0, den1p + den2p + b0, denom_eps);

                if r + EPS < best_ratio
                    best_ratio = r;
                    best_i = i;
                    best_t = t;
                end
            end
        end

        if best_i > 0
            [x, alphaA, num1, num2, alphaB, den1, den2] = apply_change_ratio( ...
                best_i, best_t, x, alphaA, num1, num2, A, a, alphaB, den1, den2, B, b);
            improved = true;
        end
    end
end

% Shake: perform k random coordinate changes (with replacement)
function [x, alphaA, num1, num2, alphaB, den1, den2] = shake_ratio( ...
    x, alphaA, num1, num2, alphaB, den1, den2, A, a, B, b, k)

    n = numel(x);
    for rep = 1:k
        i = randi(n);
        s = x(i);
        candidates = [-1,0,1];
        candidates(candidates == s) = [];
        tval = candidates(randi(2));

        [x, alphaA, num1, num2, alphaB, den1, den2] = apply_change_ratio( ...
            i, tval, x, alphaA, num1, num2, A, a, alphaB, den1, den2, B, b);
    end
end
