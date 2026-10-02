function [f_best, x_best] = vns_ratio_multistart(A, a, a0, B, b, b0, n_times, sense)

    x_best = [];
    if sense == 0
        f_best = inf;
    end
    if sense == 1
        f_best = -inf;
    end
    n = size(A, 1);

    opts.iter_max = 2;
    opts.kmin     = 2;
    opts.kmax     = n;
    opts.kstep    = 2;
    opts.verbose  = false;
    
    for i=1:n_times

        x0 = randn(n,1);
        opts.seed = i;
        if sense == 0
            [f_vns, x_vns] = vns_ternary_ratio_qp(A, a, a0, B, b, b0, x0, opts);
            if f_vns < f_best
                disp(f_vns)
                f_best = f_vns;
                x_best = x_vns;
            end
        end
        if sense == 1
            [f_vns_min, x_vns] = vns_ternary_ratio_qp(-A, -a, -a0, B, b, b0, x0, opts);
            f_vns = -f_vns_min;
            if f_vns > f_best
                disp(f_vns)
                f_best = f_vns;
                x_best = x_vns;
            end
        end
    
    end

end