function [f_best, x_best] = vns_multistart(Q, c, n_times)

    x_best = [];
    f_best = inf;
    n = size(Q, 1);

    opts.iter_max = 2;
    opts.kmin     = 2;
    opts.kmax     = n;
    opts.kstep    = 2;
    opts.verbose  = false;
    
    for i=1:n_times

        x0 = randn(n,1);
        opts.seed = i;
        [f_vns, x_vns] = vns_ternary_qp(Q, c, x0, opts);
         if f_vns < f_best
             disp(f_vns)
             f_best = f_vns;
             x_best = x_vns;
         end
    
    end

end