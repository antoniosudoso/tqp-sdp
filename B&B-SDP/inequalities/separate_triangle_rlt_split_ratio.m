function [Bcell, l, viol_vect] = separate_triangle_rlt_split_ratio(Y, eps, max_sep_ineq, max_new_ineq)

    n = size(Y, 1)-1;
    
    Bcell = cell(1, max_sep_ineq);
    l = zeros(max_sep_ineq, 1);
    viol_vect = zeros(max_sep_ineq, 1);
    c = 1;
    
    stop = false;
    for i=2:n+1
        for j=i+1:n+1

            % Xii - Xij >= 0 (TRI)
            viol = Y(i,i)-Y(i,j);
            if viol <= -eps
                %keyboard
                Bcell{c} = sparse([i, i, j], [i, j, i], [1, -0.5, -0.5], n+1, n+1);
                l(c) = 0;
                viol_vect(c) = viol;
                c = c + 1;
            end
            % Xjj - Xij >= 0 (TRI)
            viol = Y(j,j)-Y(i,j);
            if viol <= -eps
                %keyboard
                Bcell{c} = sparse([j, i, j], [j, j, i], [1, -0.5, -0.5], n+1, n+1);
                l(c) = 0;
                viol_vect(c) = viol;
                c = c + 1;
            end

            % Xii + Xij >= 0 (TRI)
            viol = Y(i,i)+Y(i,j);
            if viol <= -eps
                %keyboard
                Bcell{c} = sparse([i, i, j], [i, j, i], [1, 0.5, 0.5], n+1, n+1);
                l(c) = 0;
                viol_vect(c) = viol;
                c = c + 1;
            end
            % Xjj + Xij >= 0 (TRI)
            viol = Y(j,j)+Y(i,j);
            if viol <= -eps
                %keyboard
                Bcell{c} = sparse([j, i, j], [j, j, i], [1, 0.5, 0.5], n+1, n+1);
                l(c) = 0;
                viol_vect(c) = viol;
                c = c + 1;
            end

            
            % Xij + xi + xj >= -t (RLT)
            viol = Y(i,j)+Y(1,i)+Y(1,j)+Y(1,1);
            %keyboard
            if viol <= -eps
                %keyboard
                Bcell{c} = sparse([i, j, 1, i, 1, j, 1], [j, i, i, 1, j, 1, 1], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0], n+1, n+1);
                l(c) = 0;
                viol_vect(c) = viol;
                c = c + 1;
            end
            % Xij - xi - xj >= -t (RLT)
            viol = Y(i,j)-Y(1,i)-Y(1,j)+Y(1,1);
            if viol <= -eps
                %keyboard
                Bcell{c} = sparse([i, j, 1, i, 1, j, 1], [j, i, i, 1, j, 1, 1], [0.5, 0.5, -0.5, -0.5, -0.5, -0.5, 1.0], n+1, n+1);
                l(c) = 0;
                viol_vect(c) = viol;
                c = c + 1;
            end
            % -Xij + xi − xj >= -t (RLT)
            viol = -Y(i,j)+Y(1,i)-Y(1,j)+Y(1,1);
            if viol <= -eps
                %keyboard
                Bcell{c} = sparse([i, j, 1, i, 1, j, 1], [j, i, i, 1, j, 1, 1], [-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 1.0], n+1, n+1);
                l(c) = 0;
                viol_vect(c) = viol;
                c = c + 1;
            end
            % -Xij - xi + xj >= −t (RLT)
            viol = -Y(i,j)-Y(1,i)+Y(1,j)+Y(1,1);
            if viol <= -eps
                %keyboard
                Bcell{c} = sparse([i, j, 1, i, 1, j, 1], [j, i, i, 1, j, 1, 1], [-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 1.0], n+1, n+1);
                l(c) = 0;
                viol_vect(c) = viol;
                c = c + 1;
            end

            if c >= max_sep_ineq
                stop=true;
                break;
            end

            for k=j+1:n+1
                 if (i ~= k)

                     %  Xij + Xjk + Xik >= -t (TRI)
                     viol = Y(i,j)+Y(j,k)+Y(i,k)+Y(1,1);
                     %keyboard
                     if viol <= -eps
                         %keyboard
                         Bcell{c} = sparse([i, j, j, k, i, k, 1], [j, i, k, j, k, i, 1], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % -Xij - Xjk + Xik >= −t (TRI)
                     viol = -Y(i,j)-Y(j,k)+Y(i,k)+Y(1,1);
                     if viol <= -eps
                         %keyboard
                         Bcell{c} = sparse([i, j, j, k, i, k, 1], [j, i, k, j, k, i, 1], [-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 1.0], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % Xij - Xjk - Xik >= −t (TRI)
                     viol = Y(i,j)-Y(j,k)-Y(i,k)+Y(1,1);
                     if viol <= -eps
                         %keyboard
                         Bcell{c} = sparse([i, j, j, k, i, k, 1], [j, i, k, j, k, i, 1], [0.5, 0.5, -0.5, -0.5, -0.5, -0.5, 1.0], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % -Xij + Xjk - Xik >= -t (TRI)
                     viol = -Y(i,j)+Y(j,k)-Y(i,k)+Y(1,1);
                     if viol <= -eps
                         %keyboard
                         Bcell{c} = sparse([i, j, j, k, i, k, 1], [j, i, k, j, k, i, 1], [-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 1.0], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end

                     
                     %  Xii + Xjj + 2Xij + xi + xj ≥ 0 (SPLIT)
                     viol = Y(i,i)+Y(j,j)+2*Y(i,j)+Y(1,i)+Y(1,j);
                     if viol <= -eps
                         Bcell{c} = sparse([i, j, i, j, 1, i, 1, j], [i, j, j, i, i, 1, j, 1], ...
                             [1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % Xii + Xjj + 2Xij − xi − xj ≥ 0 (SPLIT)
                     viol = Y(i,i)+Y(j,j)+2*Y(i,j)-Y(1,i)-Y(1,j);
                     if viol <= -eps
                         Bcell{c} = sparse([i, j, i, j, 1, i, 1, j], [i, j, j, i, i, 1, j, 1], ...
                             [1, 1, 1, 1, -0.5, -0.5, -0.5, -0.5], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % Xii + Xjj − 2Xij + xi − xj ≥ 0 (SPLIT)
                     viol = Y(i,i)+Y(j,j)-2*Y(i,j)+Y(1,i)-Y(1,j);
                     if viol <= -eps
                         Bcell{c} = sparse([i, j, i, j, 1, i, 1, j], [i, j, j, i, i, 1, j, 1], ...
                             [1, 1, -1, -1, 0.5, 0.5, -0.5, -0.5], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % Xii + Xjj − 2Xij − xi + xj ≥ 0 (SPLIT)
                     viol = Y(i,i)+Y(j,j)-2*Y(i,j)-Y(1,i)+Y(1,j);
                     if viol <= -eps
                         Bcell{c} = sparse([i, j, i, j, 1, i, 1, j], [i, j, j, i, i, 1, j, 1], ...
                             [1, 1, -1, -1, -0.5, -0.5, 0.5, 0.5], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end

      
                     if c >= max_sep_ineq
                         stop=true;
                         break;
                     end
                     
                 end
            end
            
            if stop
                break;
            end
            
        end
        
        if stop
            break;
        end
        
    end
    
    n_ineq = c-1;
    if n_ineq <= max_new_ineq
        Bcell = Bcell(1:n_ineq);
        l = l(1:n_ineq);
    else    
        viol_vect = abs(viol_vect(1:n_ineq));
        [~, id_sorted] = sort(viol_vect, 'descend');
        Bcell = Bcell(id_sorted);
        l = l(id_sorted);
        Bcell = Bcell(1:max_new_ineq);
        l = l(1:max_new_ineq);
    end
    
end