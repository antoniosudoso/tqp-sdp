function result = solve_ternary_mosek_ratio(data, gub, x_gub, params)

    % This script solves the SDP relaxation of the following nonconvex qudratic programming problem:
    %
    %    min      x'*A*x + a'*x + a0 / x'*B*x + b'*x + b0
    %    s.t.     x \in {-1, 0, 1}
    %
    % * Q, c: quadratic and linear objective coefficients
    % * gub: initial global upper bound
    % * x_gub: solution attaining the initial upper bound
    % * params: structure containing the following fields
    %   - time_limit: maximum B&B running time in seconds
    %   - max_nodes: maximum number of B&B nodes, including the root
    %   - opt_gap: relative optimality-gap tolerance
    %   - n_threads: number of MOSEK threads
    %   - sdp_verbose: 0 disables solver output; 1 enables it
    %   - cp_max_iter: maximum number of iterations in each cutting-plane phase
    %   - cp_eps_ineq: minimum violation required to identify a cut
    %   - cp_max_sep_ineq: maximum number of TRI/RLT/SPLIT candidates examined
    %   - cp_max_new_tri: maximum number of TRI/RLT/SPLIT cuts added per iteration
    %   - cp_max_new_penta: maximum number of combined PENTA/EPTA cuts added
    %   - cp_ntimes_penta: number of multistarts for PENTA separation
    %   - cp_ntimes_epta: number of multistarts for EPTA separation
    
    result = struct();
    tStartBB = tic;
    tRootBB = tic;
    
    result.best_ub = gub;
    result.best_x_ub = x_gub;
    result.type = 1;

    n = data.n;

    opts.iter_max = 2;
    opts.kmin     = 2;
    opts.kmax     = n;
    opts.kstep    = 2;
    opts.verbose  = false;
    
    C_full = zeros(n+1, n+1);
    C_full(1, 1) = data.a0;
    C_full(2:n+1, 1) = 0.5 * data.a;
    C_full(1, 2:n+1) = 0.5 * data.a';
    C_full(2:n+1, 2:n+1) = data.A;
        
    X = sdpvar(n, n);
    x = sdpvar(n, 1);
    t = sdpvar(1, 1);
    X_bar = [t, x'; x, X];

    objective = C_full(:)'*X_bar(:);
    constr = [X(:)'*data.B(:) + x'*data.b + t*data.b0 == 1; ...
        X_bar >= 0; t >= 0; diag(X) >= x; diag(X) >= -x; diag(X) <= t * ones(n, 1); t <= 1];

    % solve SDP
    options_mosek = sdpsettings('verbose', params.sdp_verbose, 'solver', 'mosek', ...
        'mosek.MSK_IPAR_NUM_THREADS', double(params.n_threads));
    tStart = tic;
    optimize(constr, objective, options_mosek);
    timeSDP = toc(tStart);
    lb = value(objective);
    
    % local search
    Y = value(X_bar);
    y_SDP = Y(2:n+1, 1);  
    [f_vns, x_vns] = vns_ternary_ratio_qp(data.A, data.a, data.a0, data.B, data.b, data.b0, sign(y_SDP), opts);        
    if f_vns < result.best_ub
        result.best_ub = f_vns;
        result.best_x_ub = x_vns;
    end

    gap = (result.best_ub - lb) / max(1, abs(result.best_ub));
    result.bound_basic = lb;
    result.gap_basic = gap;
    result.best_lb = lb;

    fprintf('\n\t ********************************************\n')
    fprintf('\t # Inequalities       = %d \n', 3*n); 
    fprintf('\t Valid LB             = %.6f \n', lb);
    fprintf('\t Time SDP (s)         = %.3f \n', timeSDP);
    fprintf('\t UB                   = %.6f \n', result.best_ub);
    fprintf('\t Relative gap (%%)     = %.4f \n', gap*100);
    fprintf('\t ********************************************\n\n')

    methods = {'tri','rlt','split','penta','epta','enna'};
    for i = 1:numel(methods)
        m = methods{i};
        result.(['bound_' m]) = -inf;
        result.(['gap_' m]) = inf;
        result.(['added_' m]) = 0;
        result.(['time_' m]) = 0;
        result.(['iter_' m]) = 0;
        if ismember(m, {'penta','epta', 'enna'})
            result.(['time_sep_' m]) = 0;
        end
    end

    if result.gap_basic <= params.opt_gap, result.time_bb = toc(tRootBB); result.nodes = 1; result.gap_bb = 0; return, end

    params.cp_verbose = 1;
    [result, Y, constr] = cp_add_TRI_RLT_SPLIT_ratio(data, params, Y, X_bar, constr, objective, options_mosek, result);    
    if result.gap_tri <= params.opt_gap, result.time_bb = toc(tRootBB); result.nodes = 1; result.gap_bb = 0; return, end
    [result, Y, constr] = cp_add_PENTA_EPTA_ratio(data, params, Y, X_bar, constr, objective, options_mosek, result);
    if result.gap_penta <= params.opt_gap, result.time_bb = toc(tRootBB); result.nodes = 1; result.gap_bb = 0; return, end
    
    y_SDP = Y(2:n+1, 1);
    Y_SDP = Y(2:n+1, 2:n+1);
    t_SDP = Y(1, 1);
    %[br_idx, ~] = selectMostFractionalTernary(y_SDP ./ t_SDP);
    [br_idx, ~] = selectApproximationError(data.A, y_SDP, t_SDP .* Y_SDP);
    %approx_norm = norm(t_SDP .* Y_SDP - y_SDP*y_SDP', 'fro');
    %if br_idx == -1 && approx_norm > 1e-8
    %    [br_idx, ~] = selectApproximationError(data.A, y_SDP, t_SDP .* Y_SDP);
    %end

    queue = [];
    % explore -1, 0, 1
    newChildren = {
        struct('constr', [constr; x(br_idx) ==  t; X(br_idx, br_idx) == t], 'depth', 1, 'bound', result.best_lb);
        struct('constr', [constr; x(br_idx) ==  0; X(br_idx, br_idx) == 0], 'depth', 1, 'bound', result.best_lb);
        struct('constr', [constr; x(br_idx) == -t; X(br_idx, br_idx) == t], 'depth', 1, 'bound', result.best_lb);
    };
    queue = [queue; newChildren];

    elapsed_time = toc(tRootBB);
    nodes_processed = 0;
    global_lb = result.best_lb;

    while ~isempty(queue) && nodes_processed < (params.max_nodes - 1) && elapsed_time < params.time_limit

        elapsed_time = toc(tStartBB);

        % ---- Select node with best (lowest) relaxation bound ----
        [~, best_idx] = min(cellfun(@(node) node.bound, queue));
        node = queue{best_idx};
        queue(best_idx) = [];

        % ---- Select node with depth-first strategy (LIFO stack) ----
        %node = queue{end};
        %queue(end) = [];

        nodes_processed = nodes_processed + 1;
        fprintf(' [B&B] Processing node #%d at depth %d ...\n', nodes_processed, node.depth);

        % ---- Solve relaxation for this node ----
        tNode = tic;
        optimize(node.constr, objective, options_mosek);
        node_lb = value(objective);
        Ynode = value(X_bar);
        timeNode = toc(tNode);
        gapNode = (result.best_ub - node_lb) / max(1, abs(result.best_ub));

        params.cp_verbose = 0;

        if params.cp_verbose
            fprintf('\n\t ********************************************\n')
            fprintf('\t # Inequalities     = %d \n', length(node.constr)); 
            fprintf('\t Valid LB           = %.6f \n', node_lb);
            fprintf('\t Time SDP (s)       = %.3f \n', timeNode);
            fprintf('\t UB                 = %.6f \n', result.best_ub);
            fprintf('\t Relative gap (%%)   = %.4f \n', gapNode*100);
            fprintf('\t ********************************************\n\n')
        end

        if ~isfinite(node_lb)
            fprintf('   ✓ Pruned (infeasible) [%.3fs]\n\n', timeNode);
            continue;
        end

        % ---- Prune by bound ----
        if gapNode <= params.opt_gap
            fprintf('   ✓ Node done in %.3fs | Pruned (bound %.6f ≥ best UB %.6f)\n\n', timeNode, node_lb, result.best_ub);
            continue;
        end

        % ---- Local search & update incumbent ----
        y_relax = Ynode(2:n+1, 1);
        [f_ls, x_ls] = vns_ternary_ratio_qp(data.A, data.a, data.a0, data.B, data.b, data.b0, sign(y_relax), opts);
        if f_ls < result.best_ub
            fprintf('   * New incumbent: UB %.7f → %.7f\n', result.best_ub, f_ls);
            result.best_ub = f_ls;
            result.best_x_ub = x_ls;
            gapNode = (result.best_ub - node_lb) / max(1, abs(result.best_ub));
            if gapNode <= params.opt_gap
                fprintf('   ✓ Node done in %.3fs | Pruned (bound %.6f ≥ best UB %.6f)\n\n', timeNode, node_lb, result.best_ub);
                continue;
            end
        end

        % ------------------------- START CUTTING-PLANE -------------------------

        timeStartNodeCP = tic;
        Y0 = Ynode;
        c0 = node.constr;
        r0 = result;
        r0.best_lb = node_lb;

        [r0, Y0, c0] = cp_add_TRI_RLT_SPLIT_ratio(data, params, Y0, X_bar, c0, objective, options_mosek, r0);
        if r0.gap_tri <= params.opt_gap
            timeNode = timeNode + toc(timeStartNodeCP);
            result.best_ub = r0.best_ub;
            result.best_x_ub = r0.best_x_ub;
            fprintf('   ✓ Node done in %.3fs | Pruned (bound %.6f ≥ best UB %.6f)\n\n', timeNode, r0.best_lb, result.best_ub);
            continue;
        end
        [r0, Y0, c0] = cp_add_PENTA_EPTA_ratio(data, params, Y0, X_bar, c0, objective, options_mosek, r0);
        if r0.gap_penta <= params.opt_gap
            timeNode = timeNode + toc(timeStartNodeCP);
            result.best_ub = r0.best_ub;
            result.best_x_ub = r0.best_x_ub;
            fprintf('   ✓ Node done in %.3fs | Pruned (bound %.6f ≥ best UB %.6f)\n\n', timeNode, r0.best_lb, result.best_ub);
            continue;
        end

        node_lb = r0.best_lb;
        Ynode = Y0;
        node.constr = c0;
        timeNode = timeNode + toc(timeStartNodeCP);

        if r0.best_ub < result.best_ub
            fprintf('   * New incumbent: UB %.7f → %.7f\n', result.best_ub, r0.best_ub);
            result.best_ub = r0.best_ub;
            result.best_x_ub = r0.best_x_ub;
            gapNode = (result.best_ub - node_lb) / max(1, abs(result.best_ub));
            if gapNode <= params.opt_gap
                fprintf('   ✓ Node done in %.3fs | Pruned (bound %.6f ≥ best UB %.6f)\n\n', timeNode, node_lb, result.best_ub);
                continue;
            end
        end

        % ------------------------- END CUTTING-PLANE -------------------------

        % Update global lower bound
        global_lb = max(global_lb, node_lb);

        % ---- Branch on most-fractional variable ----
        y_relax = Ynode(2:n+1, 1);
        Y_relax = Ynode(2:n+1, 2:n+1);
        t_relax = Ynode(1, 1);
        %[br_idx, ~] = selectMostFractionalTernary(y_relax ./ t_relax);
        %approx_norm = norm(t_relax .* Y_relax - y_relax*y_relax', 'fro');
        [br_idx, ~] = selectApproximationError(data.A, y_relax, t_relax .* Y_relax);
        %if br_idx == 1 && approx_norm > 1e-8
        %    [br_idx, ~] = selectApproximationError(data.A, y_relax, t_relax .* Y_relax);
        %end
        if br_idx == -1 && approx_norm <= 1e-8
            f_int_num = x_relax'*data.A*x_relax + data.a'*x_relax + data.a0;
            f_int_den = x_relax'*data.B*x_relax + data.b'*x_relax + data.b0;
            f_int = f_int_num / f_int_den;
            if f_int < result.best_ub
                fprintf('   * Found integral feasible solution: %.6f\n', f_int);
                result.best_ub   = f_int;
                result.best_x_ub = x_relax;
            end
            fprintf('   ✓ Pruned (integral, UB=%.6f)\n\n', f_int);
            continue;
        end
        fprintf('   ↳ Branching on variable %d\n', br_idx);

        % Generate children (3-way split)
        newChildren = {
            struct('constr', [node.constr; x(br_idx) ==   t; X(br_idx, br_idx) == t], 'depth', node.depth+1, 'bound', node_lb);
            struct('constr', [node.constr; x(br_idx) ==   0; X(br_idx, br_idx) == 0], 'depth', node.depth+1, 'bound', node_lb);
            struct('constr', [node.constr; x(br_idx) ==  -t; X(br_idx, br_idx) == t], 'depth', node.depth+1, 'bound', node_lb);
        };
        queue = [queue; newChildren];

        % ---- Progress log ----
        fprintf('   ✓ Node done in %.3fs | Current LB=%.6f | Best LB=%.6f | GUB=%.6f | Gap=%.4f%% | Queue=%d\n\n', ...
            timeNode, node_lb, global_lb, result.best_ub, 100*(result.best_ub - global_lb)/max(1, abs(result.best_ub)), length(queue));

    end

    % ---- Final results ----
    if isempty(queue)
        result.gap_bb = 0;
    else
        [min_lb, ~] = min(cellfun(@(node) node.bound, queue));
        result.gap_bb   = (result.best_ub - min_lb) / max(1, abs(result.best_ub));
    end
    result.nodes    = nodes_processed + 1;
    result.time_bb  = toc(tStartBB);


end

