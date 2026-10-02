function [result, Y, constr] = cp_add_TRI_RLT_SPLIT_ratio(data, params, Y, X_bar, constr, objective, options_mosek, result_old)

    n = data.n;

    result = result_old;
    result.bound_tri = -inf;
    result.gap_tri = inf;
    result.added_tri = 0;
    result.time_tri = 0;
    result.iter_tri = 0;
    tStartTot = tic;

    for i=1:params.cp_max_iter
        
        [B_triangle, l_triangle, ~] = separate_triangle_rlt_split_ratio(Y, params.cp_eps_ineq, params.cp_max_sep_ineq, params.cp_max_new_tri);
        n_triangle = size(B_triangle, 2);
        if n_triangle <= round(n/2)
            break;
        end

        result.added_tri = result.added_tri + n_triangle;
        result.iter_tri = result.iter_tri + 1;
       
        for j=1:n_triangle
            constr = [constr; B_triangle{j}(:)'*X_bar(:) >= l_triangle(j)];
        end
        
        % solve SDP
        tStart = tic;
        optimize(constr, objective, options_mosek);
        timeSDP = toc(tStart);      
        lb = value(objective);

        Y = value(X_bar);
    
        result.bound_tri = lb;
        result.gap_tri = (result.best_ub - lb) / max(1, abs(result.best_ub));
        result.best_lb = lb;

        if params.cp_verbose
            fprintf('\n\t ********************************************\n')
            fprintf('\t Iter                            = %d \n', i); 
            fprintf('\t # TRI+RLT+SPLIT Inequalities    = %d \n', size(B_triangle, 2)); 
            fprintf('\t Valid LB                        = %.6f \n', lb);
            fprintf('\t Time SDP (s)                    = %.3f \n', timeSDP);
            fprintf('\t UB                              = %.6f \n', result.best_ub);
            fprintf('\t Relative gap (%%)                = %.4f \n', result.gap_tri*100);
            fprintf('\t ********************************************\n\n')
        end

        if result.gap_tri <= params.opt_gap
            break;
        end

    end

    result.time_tri = toc(tStartTot);



end