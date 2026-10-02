function [idx, value] = selectApproximationError(Q, x, X)

    tol = 1e-8;

    n = size(X, 1);
    values = zeros(n, 1);
    for i=1:n
        sum = 0;
        for j=1:n
            sum = sum + abs(Q(i,j)*(x(j)*x(i)-X(i,j)));
        end
        values(i) = sum;
    end
    
    pairs = sortrows([(1:n)', values], 2, 'descend'); % pairs (id_x, viol)
    %disp(pairs)
    value = pairs(1, 2);

    if value <= tol
        idx = -1;
    else
        idx = pairs(1, 1);
    end

end


