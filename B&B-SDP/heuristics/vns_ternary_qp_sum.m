function [best_cost, x_best] = vns_ternary_qp_sum(A, a, x0, opts)
% VNS_TERNARY_QP_SUM
%   Variable Neighborhood Search for TQP-Linear
%       min  x' * A * x + a' * x
%       s.t. sum(x) = 0
%            x_i in {-1, 0, 1}
%
%   The local-search neighborhood and the shaking phase use feasible paired
%   changes on distinct indices i and j. If delta is applied to x(i), then
%   -delta is applied to x(j), so sum(x) remains zero after every move.
%
% INPUT
%   A    : (n x n) symmetric matrix
%   a    : (n x 1) vector
%   x0   : (n x 1) initial vector (projected to a balanced ternary vector)
%   opts : struct with optional fields:
%            .iter_max        (default 3)
%            .kmin            (default 2)
%            .kmax            (default n)
%            .kstep           (default 2)
%            .seed            (default 1)
%            .max_local_iters (default 1e2)
%            .verbose         (default false)
%
% OUTPUT
%   x_best   : best solution found, with entries in {-1,0,1} and sum zero
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

    % ----- Options with the values used in the paper -----
    iter_max        = get_opt(opts, 'iter_max',        3);
    kmin            = get_opt(opts, 'kmin',            2);
    kmax            = get_opt(opts, 'kmax',            n);
    kstep           = get_opt(opts, 'kstep',           2);
    seed            = get_opt(opts, 'seed',            1);
    max_local_iters = get_opt(opts, 'max_local_iters', 1e3);
    verbose         = get_opt(opts, 'verbose',         false);

    if iter_max < 0 || iter_max ~= floor(iter_max)
        error('opts.iter_max must be a nonnegative integer.');
    end
    if kmin < 1 || kstep < 1 || kmin ~= floor(kmin) || kstep ~= floor(kstep)
        error('opts.kmin and opts.kstep must be positive integers.');
    end
    if kmax < 0 || kmax ~= floor(kmax)
        error('opts.kmax must be a nonnegative integer.');
    end
    if max_local_iters < 0 || max_local_iters ~= floor(max_local_iters)
        error('opts.max_local_iters must be a nonnegative integer.');
    end

    % ----- Initial feasible solution -----
    % For each possible common number of +1 and -1 entries, choose the
    % largest components of x0 as +1 and the smallest as -1. Select the
    % closest such balanced ternary vector in Euclidean distance.
    x = project_balanced_ternary(double(x0(:)));

    % ----- Initialize aggregates -----
    alphaA = A * x;
    f1     = x' * A * x;
    f2     = a' * x;

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

            % ---- Shake: random feasible paired changes ----
            [x, f1, f2, alphaA] = shake_sum(x, f1, f2, alphaA, A, a, k);

            % ---- Local search: best feasible paired change ----
            [x, f1, f2, alphaA] = local_search_sum( ...
                x, f1, f2, alphaA, A, a, max_local_iters);

            cur_cost = f1 + f2;

            if cur_cost + EPS < best_cost
                % Improvement: update best.
                x_best    = x;
                alphaA_b  = alphaA;
                f1_b      = f1;
                f2_b      = f2;
                best_cost = cur_cost;

                if verbose
                    fprintf('Improved: it=%d k=%d obj=%.6g\n', it, k, best_cost);
                end

                % Intensify.
                k = kmin;
            else
                % No improvement: restore best and diversify.
                x      = x_best;
                alphaA = alphaA_b;
                f1     = f1_b;
                f2     = f2_b;
                k      = k + kstep;
            end
        end
    end

    % Remove any accumulated floating-point error in the reported value.
    best_cost = x_best' * A * x_best + a' * x_best;
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

% Euclidean projection onto {x in {-1,0,1}^n : sum(x) = 0}.
function x = project_balanced_ternary(y)
    n = numel(y);
    [y_sorted, order] = sort(y, 'ascend');

    prefix = [0; cumsum(y_sorted)];
    baseline = y' * y;
    best_projection_cost = baseline;
    best_m = 0;

    % A balanced ternary vector has the same number m of +1 and -1 entries.
    for m = 1:floor(n / 2)
        sum_smallest = prefix(m + 1);
        sum_largest  = prefix(n + 1) - prefix(n - m + 1);
        projection_cost = baseline + 2 * m ...
            + 2 * sum_smallest - 2 * sum_largest;

        if projection_cost < best_projection_cost
            best_projection_cost = projection_cost;
            best_m = m;
        end
    end

    x = zeros(n, 1);
    if best_m > 0
        x(order(1:best_m)) = -1;
        x(order(n-best_m+1:n)) = 1;
    end
end

% Objective after the feasible paired move
%   x(i) := x(i) + delta, x(j) := x(j) - delta.
function [f1p, f2p] = trial_pair_change( ...
        i, j, delta, f1, f2, alphaA, A, a)

    delta_i = delta;
    delta_j = -delta;

    f1p = f1 ...
        + 2 * delta_i * alphaA(i) ...
        + 2 * delta_j * alphaA(j) ...
        + delta_i^2 * A(i,i) ...
        + delta_j^2 * A(j,j) ...
        + 2 * delta_i * delta_j * A(i,j);
    f2p = f2 + delta_i * a(i) + delta_j * a(j);
end

% Apply a feasible paired move and update A*x.
function [x, f1, f2, alphaA] = apply_pair_change( ...
        i, j, delta, x, f1, f2, alphaA, A, a)

    [f1, f2] = trial_pair_change(i, j, delta, f1, f2, alphaA, A, a);

    delta_i = delta;
    delta_j = -delta;
    alphaA = alphaA + delta_i * A(:,i) + delta_j * A(:,j);
    x(i) = x(i) + delta_i;
    x(j) = x(j) + delta_j;
end

% Best-improvement local search in the feasible paired-move neighborhood.
function [x, f1, f2, alphaA] = local_search_sum( ...
        x, f1, f2, alphaA, A, a, max_iters)

    n = numel(x);
    EPS = 1e-6;
    it_count = 0;
    deltas = [-2, -1, 1, 2];

    improved = true;
    while improved && it_count < max_iters
        it_count = it_count + 1;

        improved  = false;
        best_val  = f1 + f2;
        best_i    = 0;
        best_j    = 0;
        best_delta = 0;

        for i = 1:n-1
            for j = i+1:n
                for d_idx = 1:numel(deltas)
                    delta = deltas(d_idx);
                    value_i = x(i) + delta;
                    value_j = x(j) - delta;

                    if value_i < -1 || value_i > 1 || ...
                            value_j < -1 || value_j > 1
                        continue;
                    end

                    [f1p, f2p] = trial_pair_change( ...
                        i, j, delta, f1, f2, alphaA, A, a);
                    val = f1p + f2p;
                    if val + EPS < best_val
                        best_val   = val;
                        best_i     = i;
                        best_j     = j;
                        best_delta = delta;
                    end
                end
            end
        end

        if best_i > 0
            [x, f1, f2, alphaA] = apply_pair_change( ...
                best_i, best_j, best_delta, x, f1, f2, alphaA, A, a);
            improved = true;
        end
    end
end

% Shake by applying random feasible paired moves. The radius k counts the
% nominal number of perturbed coordinates, so ceil(k/2) pairs are applied.
function [x, f1, f2, alphaA] = shake_sum( ...
        x, f1, f2, alphaA, A, a, k)

    n = numel(x);
    if n < 2
        return;
    end

    n_pair_moves = ceil(k / 2);
    deltas = [-2, -1, 1, 2];

    for move = 1:n_pair_moves
        found = false;

        % Rejection sampling normally finds a feasible pair immediately.
        for attempt = 1:max(20, 2 * n)
            pair = randperm(n, 2);
            i = pair(1);
            j = pair(2);
            feasible = deltas( ...
                x(i) + deltas >= -1 & x(i) + deltas <= 1 & ...
                x(j) - deltas >= -1 & x(j) - deltas <= 1);

            if ~isempty(feasible)
                delta = feasible(randi(numel(feasible)));
                found = true;
                break;
            end
        end

        % Deterministic fallback for the unlikely case of repeated rejection.
        if ~found
            for i = 1:n-1
                for j = i+1:n
                    feasible = deltas( ...
                        x(i) + deltas >= -1 & x(i) + deltas <= 1 & ...
                        x(j) - deltas >= -1 & x(j) - deltas <= 1);
                    if ~isempty(feasible)
                        delta = feasible(randi(numel(feasible)));
                        found = true;
                        break;
                    end
                end
                if found
                    break;
                end
            end
        end

        if ~found
            return;
        end

        [x, f1, f2, alphaA] = apply_pair_change( ...
            i, j, delta, x, f1, f2, alphaA, A, a);
    end
end
