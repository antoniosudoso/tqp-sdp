function [best_cost, x_best] = vns_ternary_qp(A, a, x0, opts)
% VNS_TERNARY_QP
%   Variable Neighborhood Search for
%       min  x' * A * x + a' * x
%       s.t. x_i in {-1, 0, 1}
%
% INPUT
%   A    : (n x n) symmetric matrix
%   a    : (n x 1) vector
%   x0   : (n x 1) initial vector (projected to {-1,0,1})
%   opts : struct with optional fields:
%            .iter_max        (default 3)
%            .kmin            (default 2)
%            .kmax            (default n)
%            .kstep           (default 2)
%            .seed            (default 1)
%            .max_local_iters (default 1e2)  % safety against very long LS
%            .verbose         (default false)
%
% OUTPUT
%   x_best   : best solution found (n x 1, entries in {-1,0,1})
%   best_cost: objective value at x_best


    if nargin < 4
        opts = struct();
    end

    % ----- Basic checks -----
    [n1, n2] = size(A);
    if n1 ~= n2
        error('A must be square.');
    end
    n = n1;
    if numel(a) ~= n
        error('Dimension mismatch: a must have length size(A,1).');
    end
    if numel(x0) ~= n
        error('Dimension mismatch: x0 must have length size(A,1).');
    end

    A = double(A);
    a = double(a(:));

    % ----- Options with defaults -----
    iter_max        = get_opt(opts, 'iter_max',        3);
    kmin            = get_opt(opts, 'kmin',            2);
    kmax            = get_opt(opts, 'kmax',            n);
    kstep           = get_opt(opts, 'kstep',           2);
    seed            = get_opt(opts, 'seed',            1);
    max_local_iters = get_opt(opts, 'max_local_iters', 1e3);
    verbose         = get_opt(opts, 'verbose',         false);

    % ----- Initial solution and projection to {-1,0,1} -----
    x = double(x0(:));
    x = round(x);
    x(x >  1) =  1;
    x(x < -1) = -1;

    % ----- Initialize aggregates -----
    alphaA = A * x;      % A*x
    f1     = x' * A * x; % quadratic term
    f2     = a' * x;     % linear term

    x_best    = x;
    alphaA_b  = alphaA;
    f1_b      = f1;
    f2_b      = f2;
    best_cost = f1 + f2;

    % ----- RNG -----
    rng(seed);

    EPS = 1e-6;
    % ===================== MAIN VNS LOOP =====================
    for it = 1:iter_max
        k = kmin;
        while k <= kmax

            % ---- Shake: random changes on k coordinates ----
            [x, f1, f2, alphaA] = shake(x, f1, f2, alphaA, A, a, k);

            % ---- Local search: best 1-change moves ----
            [x, f1, f2, alphaA] = local_search(x, f1, f2, alphaA, A, a, max_local_iters);

            cur_cost = f1 + f2;

            if cur_cost + EPS < best_cost
                % Improvement: update best
                x_best    = x;
                alphaA_b  = alphaA;
                f1_b      = f1;
                f2_b      = f2;
                best_cost = cur_cost;

                if verbose
                    fprintf('Improved: it=%d k=%d obj=%.6g\n', it, k, best_cost);
                end

                % Intensify
                k = kmin;
            else
                % No improvement: restore best & diversify
                x      = x_best;
                alphaA = alphaA_b;
                f1     = f1_b;
                f2     = f2_b;
                k      = k + kstep;
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

% Current objective
function obj = cur_obj(f1, f2)
    obj = f1 + f2;
end

% Trial change x_i := t (t in {-1,0,1}, t ~= x_i)
% Correct formula:
%   delta = t - s
%   f1' = f1 + 2*delta*alphaA(i) + delta^2 * A(i,i)
%   f2' = f2 + delta * a(i)
function [f1p, f2p] = trial_change(i, t, x, f1, f2, alphaA, A, a)
    s = x(i);
    if t == s
        f1p = f1;
        f2p = f2;
        return;
    end
    delta = t - s;
    Aii   = A(i,i);

    f1p = f1 + 2 * delta * alphaA(i) + (delta^2) * Aii;
    f2p = f2 + delta * a(i);
end

% Apply the change x_i := t (update x, f1, f2, alphaA)
function [x, f1, f2, alphaA] = apply_change(i, t, x, f1, f2, alphaA, A, a)
    s = x(i);
    if t == s
        return;
    end
    delta = t - s;
    Aii   = A(i,i);

    % Update objective components
    f1 = f1 + 2 * delta * alphaA(i) + (delta^2) * Aii;
    f2 = f2 + delta * a(i);

    % Update alphaA = A*x
    alphaA = alphaA + delta * A(:,i);

    % Update x
    x(i) = t;
end

% Best-improvement local search in 1-change neighborhood
function [x, f1, f2, alphaA] = local_search(x, f1, f2, alphaA, A, a, max_iters)
    n   = numel(x);
    EPS = 1e-6;
    it_count = 0;

    improved = true;
    while improved && it_count < max_iters
        it_count = it_count + 1;

        improved = false;
        best_val = cur_obj(f1, f2);
        best_i   = 0;
        best_t   = 0;

        for i = 1:n
            s = x(i);
            % the two alternative values in {-1,0,1} \ {s}
            candidates = [-1, 0, 1];
            candidates(candidates == s) = [];

            for c_idx = 1:2
                t = candidates(c_idx);
                [f1p, f2p] = trial_change(i, t, x, f1, f2, alphaA, A, a);
                val = f1p + f2p;
                if val + EPS < best_val
                    best_val = val;
                    best_i   = i;
                    best_t   = t;
                end
            end
        end

        if best_i > 0
            [x, f1, f2, alphaA] = apply_change(best_i, best_t, ...
                                               x, f1, f2, alphaA, A, a);
            improved = true;
        end
    end
end

% Shake: perform k random coordinate changes (with replacement)
function [x, f1, f2, alphaA] = shake(x, f1, f2, alphaA, A, a, k)
    n = numel(x);
    for t = 1:k
        i = randi(n);
        s = x(i);
        candidates = [-1, 0, 1];
        candidates(candidates == s) = [];
        tval = candidates(randi(2));
        [x, f1, f2, alphaA] = apply_change(i, tval, x, f1, f2, alphaA, A, a);
    end
end
