function [result, Y, constr] = cp_add_PENTA_EPTA_ratio(data, params, Y, X_bar, constr, objective, options_mosek, result_old)

    n = data.n;

    result = result_old;
    result.bound_penta = -inf;
    result.gap_penta = inf;
    result.added_penta = 0;
    result.time_penta = 0;
    result.iter_penta = 0;
    result.time_sep_penta = 0;
    tStartTot = tic;

    for i=1:params.cp_max_iter
        tStartSep = tic;
        [B_penta, l_penta, ~] = separate_penta_epta_multi_ratio(Y, params.cp_eps_ineq, params.cp_max_new_penta, params.cp_ntimes_penta, params.cp_ntimes_epta);
        result.time_sep_penta = result.time_sep_penta + toc(tStartSep);
        n_penta = size(B_penta, 2);
        if n_penta <= n
            break;
        end

        result.added_penta = result.added_penta + n_penta;
        result.iter_penta = result.iter_penta + 1;

        for j=1:n_penta
            constr = [constr; B_penta{j}(:)'*X_bar(:) >= l_penta(j)];
        end
        
        % solve SDP
        tStart = tic;
        optimize(constr, objective, options_mosek);
        timeSDP = toc(tStart);      
        lb = value(objective);
        
        Y = value(X_bar);

        result.bound_penta = lb;
        result.gap_penta = (result.best_ub - lb) / max(1, abs(result.best_ub));
        result.best_lb = lb;

        if params.cp_verbose
            fprintf('\n\t ********************************************\n')
            fprintf('\t Iter                         = %d \n', i); 
            fprintf('\t # PENTA+EPTA Inequalities    = %d \n', size(B_penta, 2)); 
            fprintf('\t Valid LB                     = %.6f \n', lb);
            fprintf('\t Time SDP (s)                 = %.3f \n', timeSDP);
            fprintf('\t UB                           = %.6f \n', result.best_ub);
            fprintf('\t Relative gap (%%)             = %.4f \n', result.gap_penta*100);
            fprintf('\t ********************************************\n\n')
        end

        if result.gap_penta <= params.opt_gap
            break;
        end

    end

     result.time_penta = toc(tStartTot);


end