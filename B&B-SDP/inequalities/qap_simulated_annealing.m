function [cost, ineq, perm_best] = qap_simulated_annealing(H, k, X)
    
    %QAP_SIMULATED_ANNEALING  Simulated annealing heuristic for:
    %   min_P  < (P X P^T)(1:k,1:k), H >  over permutation matrices P
    %
    % Inputs
    %   H  : k x k matrix (defines the k-gonal inequality)
    %   k  : size of the submatrix (e.g., 5 for pentagonal, 7 for heptagonal)
    %   X  : n x n real matrix (will be symmetrized by copying upper -> lower)
    %
    % Outputs
    %   cost      : best objective value found
    %   ineq      : 1-by-k integer vector with the first k entries of the best
    %               permutation
    %   perm_best : 1-by-n best full permutation found
    %
    % Notes
    % - This is a direct MATLAB implementation of the provided C code by T. Hrga. The acceptance rule,
    %   temperature reduction (0.6), and inner-iteration growth (x1.1) are unchanged.

    % Basic checks
    if size(H,1) ~= k || size(H,2) ~= k
        error('H must be k x k.');
    end
    [n1, n2] = size(X);
    if n1 ~= n2
        error('X must be square.');
    end
    n = n1;

    % Copy upper triangle to lower (make symmetric, keep upper as-is)
    % X = triu(X) + triu(X,1)';

    % Parameters
    inner_iter    = n;
    reduce_temp   = 0.6;
    increase_iter = 1.1;

    % Initial temperature:
    % t = (sum(sum(abs(H))) * sum(sum(abs(X)))) / (n*(n-1))
    % The C code replaces sum(abs(H)) by k^2.
    t = k^2;
    sumX = trace(X) + 2*sum(sum(abs(triu(X,1))));
    t = t * sumX / (n*(n-1));

    % Best found cost
    cost = inf;

    % Initial random permutation (MATLAB 1..n)
    perm = randperm(n);
    perm_best = perm;

    % Initialize temperature and inner-iterations
    t1 = t;
    m1 = inner_iter;

    % Compute current solution sol = <H, (P X P^T)(1:k,1:k)>
    % Implemented as in the C code:
    %   sol = sum_i X(perm(i),perm(i)) + sum_{i<j} 2*H(j,i)*X(perm(j),perm(i))
    sol = 0.0;
    for i = 1:k
        sol = sol + X(perm(i), perm(i));
        for j = i+1:k
            sol = sol + 2 * H(j,i) * X(perm(j), perm(i));
        end
    end

    % Main annealing loop: continue until a full outer pass makes no real change
    while true
        not_done = false;

        % m1 iterations at constant temperature
        for num_it = 1:m1
            % Pick i1 in {1..k} and i2 in {1..n}, then enforce i1 <= i2
            i1 = randi(k);
            i2 = randi(n);
            if i2 < i1
                tmp = i1; i1 = i2; i2 = tmp;
            end

            % Delta from swapping perm(i1) and perm(i2)
            delta = 0.0;

            % First part: rows/cols interacting with position i1
            for i = 1:k
                delta = delta + H(i, i1) * ( X(perm(i), perm(i2)) - X(perm(i), perm(i1)) );
            end

            % If i2 is within the first k block, also account for its interactions
            if i2 <= k
                for i = 1:k
                    delta = delta + H(i, i2) * ( X(perm(i), perm(i1)) - X(perm(i), perm(i2)) );
                end
            end

            delta = 2.0 * delta;

            % Diagonal/self-term adjustment
            tempValue1 = X(perm(i1), perm(i1)) + X(perm(i2), perm(i2)) - 2.0 * X(perm(i2), perm(i1));

            tempValue2 = 1.0;
            if i2 <= k
                tempValue2 = tempValue2 + (-2 * H(i1, i2) + 1.0);
            end

            delta = delta + tempValue1 * tempValue2;

            sol_temp = sol + delta;

            % Acceptance rule
            if delta > 0
                dt1 = delta / t1;
                if dt1 > 5
                    accept = false;
                else
                    prob = exp(-dt1);
                    accept = (rand < prob);
                end
            else
                accept = true;
            end

            % Apply swap if accepted
            if accept
                if abs(delta) > 1e-4
                    not_done = true;
                end

                % swap perm(i1) and perm(i2)
                tmp = perm(i1); perm(i1) = perm(i2); perm(i2) = tmp;

                sol = sol_temp;

                % Track best
                if sol < cost
                    %disp(cost)
                    cost = sol;
                    perm_best = perm;
                end
            end
        end

        % Anneal
        t1 = t1 * reduce_temp;
        m1 = round(m1 * increase_iter);

        % Stop if nothing meaningful happened in this outer loop
        if ~not_done
            break;
        end
    end

    % Return first k entries of best permutation as the k-gonal inequality
    ineq = perm_best(1:k);
end
