function [Bcell, l, viol_vect] = separate_penta_epta_multi_ratio(Y, eps, max_ineq, n_times_penta, n_times_epta)

    Bcell = [];
    l = [];
    viol_vect = [];
    
    [BcellPENTA, lPENTA, viol_vectPENTA] = separate_kgonal_ratio(Y, 5, eps, max_ineq, n_times_penta);
    Bcell = [Bcell, BcellPENTA];
    l = [l; lPENTA];
    viol_vect = [viol_vect; viol_vectPENTA];

    [BcellEPTA, lEPTA, viol_vectEPTA] = separate_kgonal_ratio(Y, 7, eps, max_ineq, n_times_epta);
    Bcell = [Bcell, BcellEPTA];
    l = [l; lEPTA];
    viol_vect = [viol_vect; viol_vectEPTA];

    n_ineq = length(viol_vect);
    if n_ineq <= max_ineq
        Bcell = Bcell(1:n_ineq);
        l = l(1:n_ineq);
        viol_vect = viol_vect(1:n_ineq);
    else
        viol_vect = abs(viol_vect(1:n_ineq));
        [~, id_sorted] = sort(viol_vect, 'descend');
        %figure(1);
        %plot(viol_sorted);
        Bcell = Bcell(id_sorted);
        l = l(id_sorted);
        Bcell = Bcell(1:max_ineq);
        l = l(1:max_ineq);
    end


    
end