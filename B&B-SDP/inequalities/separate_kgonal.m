function [Bcell, l, viol_vect] = separate_kgonal(Y, k, eps, max_sep_ineq, n_times)


    n = size(Y, 1)-1;
    
    Bcell = cell(1, max_sep_ineq);
    l = zeros(max_sep_ineq, 1);
    viol_vect = zeros(max_sep_ineq, 1);
    c = 1;

    if k == 5
        rhs = 2;
        e1 = [1; 1; 1; 1; 1];
        e2 = [-1; 1; 1; 1; 1];
        e3 = [-1; -1; 1; 1; 1];
        eCell = cell(1, 3);
        eCell{1} = e1;
        eCell{2} = e2;
        eCell{3} = e3;
        nH = 3;
    elseif k == 7
        rhs = 3;
        e1 = [1; 1; 1; 1; 1; 1; 1];
        e2 = [-1; 1; 1; 1; 1; 1; 1];
        e3 = [-1; -1; 1; 1; 1; 1; 1];
        e4 = [-1; -1; -1; 1; 1; 1; 1];
        eCell = cell(1, 4);
        eCell{1} = e1;
        eCell{2} = e2;
        eCell{3} = e3;
        eCell{4} = e4;
        nH = 4;
    elseif k == 9
        rhs = 4;
        e1 = [1; 1; 1; 1; 1; 1; 1; 1; 1];
        e2 = [-1; 1; 1; 1; 1; 1; 1; 1; 1];
        e3 = [-1; -1; 1; 1; 1; 1; 1; 1; 1];
        e4 = [-1; -1; -1; 1; 1; 1; 1; 1; 1];
        e5 = [-1; -1; -1; -1; 1; 1; 1; 1; 1];
        eCell = cell(1, 5);
        eCell{1} = e1;
        eCell{2} = e2;
        eCell{3} = e3;
        eCell{4} = e4;
        eCell{5} = e5;
        nH = 5;
    else
        error('\n Invalid value for k: %d. Expected 5 - 7 - 9\n', k);
    end

    uniqueVectors = [];

    X = Y(2:n+1, 2:n+1);

    for s=1:n_times

        for t=1:nH
    
            e = eCell{t};
            H = e*e';
            [~, ineq, ~] = qap_simulated_annealing(H, k, X);
            lhs = 0;
            I = [];
            J = [];
            V = [];
            for i=1:k
                for j=i+1:k
                    lhs = lhs + (e(i)*e(j)*X(ineq(i), ineq(j)));
                    I = [I; ineq(i)+1; ineq(j)+1];
                    J = [J; ineq(j)+1; ineq(i)+1];
                    V = [V; 0.5*e(i)*e(j); 0.5*e(j)*e(i)];
                end
            end
            viol = lhs + rhs;
            %keyboard
            ineq_sorted = sort(ineq);
            if viol <= -eps && (isempty(uniqueVectors) || ~ismember(ineq_sorted, uniqueVectors, 'rows'))
                uniqueVectors = [uniqueVectors; ineq_sorted];
                Bcell{c} = sparse(I, J, V, n+1, n+1);
                l(c) = -rhs;
                viol_vect(c) = viol;
                c = c + 1;
            end
    
        end

        if c >= max_sep_ineq
            break;
        end

    end

    n_ineq = c - 1;
    Bcell = Bcell(1:n_ineq);
    l = l(1:n_ineq);
    viol_vect = viol_vect(1:n_ineq);

    
end