function result = compute_dinkelbach_our(lambda_init, data, tol, params, sense)
    
    % Dinkelbach's algorithm for fractional quadratic objective.
    %
    %
    % Inputs
    %   lambda_init : initial lambda (scalar)
    %   data        : struct with fields:
    %                 n (int)
    %                 A (n x n), a (n x 1), a0 (scalar)
    %                 B (n x n), b (n x 1), b0 (scalar)
    %   tol         : stopping tolerance on subproblem optimum value (scalar)
    %
    % Output (struct)
    %   result.x_opt        : optimal x (n x 1)
    %   result.lambda       : final lambda (scalar)
    %   result.iter         : iterations used (int)
    %   result.elapsed_time : elapsed time in seconds (double)

    fprintf('\t Running Dinkelbach''s algorithm...\n');

    n = data.n;
    x_opt = zeros(n, 1);

    iter = 0;
    lambda = lambda_init;

    tStart = tic;

    result = struct();
    result.node_list = [];
    result.time_list = [];

    while true

        iter = iter + 1;

        % ---- Solve subproblem ----
        if sense == 0
            Q = data.A - lambda .* data.B;
            c = data.a - lambda .* data.b; 
            const = data.a0 - lambda .* data.b0;
        else
            Q = lambda .* data.B - data.A;
            c = lambda .* data.b - data.a;
            const = lambda .* data.b0 - data.a0;
        end

        %keyboard

        disp('RUNNING INNER VNS...')
        [best_ub, best_x_ub] = vns_multistart(Q, c, 100);
        disp("INNER VNS DONE")

        r = solve_ternary_mosek_unc(Q, c, best_ub, best_x_ub, params);
        result.time_list = [result.time_list; r.time_bb];
        result.node_list = [result.node_list; r.nodes];
        %keyboard

        if r.time_bb > params.time_limit
            break
        end

        obj_min = r.best_ub;
        x_opt = r.best_x_ub;
        obj_min = obj_min + const;

        if sense == 0
            phi = obj_min;
        else
            phi = -obj_min;
        end

        fprintf('\t ITER %d\t\tPHI: %.12g\t\tLAMBDA: %.12g\n', iter, phi, lambda);

        if abs(phi) < tol
            break;
        end

        % Compute numerator / denominator at x_opt
        num = x_opt'*(data.A*x_opt) + x_opt'*data.a + data.a0;
        den = x_opt'*(data.B*x_opt) + x_opt'*data.b + data.b0;

        % Basic safety check
        if den == 0
            error('Denominator evaluated to 0 at iter %d; cannot update lambda.', iter);
        end

        % Update lambda
        lambda = num / den;

    end

    elapsed_time = toc(tStart);

    fprintf('\t Optimal solution value: %.12g\n', lambda);

    result.x_opt = x_opt;
    result.lambda = lambda;
    result.iter = iter;
    result.elapsed_time = elapsed_time;

end
