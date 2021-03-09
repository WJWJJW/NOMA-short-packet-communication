% clc; clear variables; close all;
% c_matrix = [600 670 960 560
%      900 280 970 540
%      310 350 950 820
%      325 290 600 540];
% c_matrix = [2 15 13 4;10 4 14 15;9 14 16 13;7 8 11 9]
% c_matrix = [108 125 150;150 135 175;122 148 250]
% c_matrix = [1 2 3;2 4 6;3 6 9]

function [starZ] = Hungarian_algorithm(c_matrix,K)       
      
    % Step 1: Row reduction
    % Row reduction
   
    row_min = min(c_matrix,[],2);
    c_matrix_p = c_matrix - row_min;
    
    
    % Step 2: Find index of element = 0
    z_matrix = (c_matrix_p == 0);

    % row scaning, find the first zero (star zero) element in each row, and ignore other 0
    % in that cloumn
    % store the column index of star zero in the starZ

    starZ = zeros(K);
    for ii = 1:K
        if find(z_matrix(ii,:)) ~= 0
            c = find(z_matrix(ii,:));
            starZ(ii, c(1)) = true;
            z_matrix(ii,:)=false;
            z_matrix(:,c)=false;
        end
    end


    primeZ=zeros(K);
    cover_col = false(1,K);
    cover_row = false(K,1);
    while 1
        if all(starZ>0)
            break;
        end

        % Step 3: Cover each column with star zero
        cover_col(sum(starZ) == 1) = true;
%         check_row(sum(starZ') == 1) = true;

        % if # of cover column = K, solution find
        if sum(cover_col) == K
            break;
        end

        % Step 4: Adjust the solution.
        while 1
            % find zero at non cover row or column, and prime it.
            if any(c_matrix_p(~cover_row,~cover_col)==0,'all') == true
                [r_prime_z,c_prime_z] = find(c_matrix_p(~cover_row,~cover_col)==0,1);
                % index tranformation
                c_row = find(~cover_row);
                c_col = find(~cover_col);
                r_prime_z = c_row(r_prime_z);
                c_prime_z = c_col(c_prime_z);
                primeZ(r_prime_z,c_prime_z) = 1;
                % check there are star zero at the row of prime zero
                if any(starZ(r_prime_z,:)) % Yes
                    % uncover column
                    cover_col(starZ(r_prime_z,:) == 1) = false;
                    % cover row
                    cover_row(r_prime_z) = true;
                else % No
                    % Step 5: Update star zero
%                     t = 0;
                    while 1
%                         if t == 0
%                             [~,c_prime_z] = find(primeZ==1,1);
%                             [r_prime_z,c_prime_z] = find(primeZ(:,c_prime_z)==1);
%                             r_prime_z = max(r_prime_z);
%                             c_prime_z = c_prime_z(1);
%                         else
%                             [~,c_prime_z] = find(primeZ(r_star_z,:)==1,1);
%                             r_prime_z = r_star_z;
%                         end
%                         primeZ(r_prime_z,c_prime_z) = 0;
%                         if any(starZ(:,c_prime_z)) == 1
%                             [r_star_z,~] = find(starZ(:,c_prime_z)==1,1);
%                             starZ(r_star_z,c_prime_z) = 0;
%                             starZ(r_prime_z,c_prime_z) = 1;
%                         else
%                             starZ(r_prime_z,c_prime_z) = 1;
%                             cover_col = false(1,K);
%                             cover_row = false(K,1);
%                             primeZ=zeros(K);
%                             break;
%                         end
%                         t = t+1;
                        if any(starZ(:,c_prime_z)) == 1
                            [r_star_z,~] = find(starZ(:,c_prime_z)==1,1);
                            starZ(r_star_z,c_prime_z) = 0;
                            starZ(r_prime_z,c_prime_z) = 1;
                        else
                            starZ(r_prime_z,c_prime_z) = 1;
                            cover_col = false(1,K);
                            cover_row = false(K,1);
                            primeZ=zeros(K);
                            break;
                        end
                        [~,c_prime_z] = find(primeZ(r_star_z,:)==1,1);
                        r_prime_z = r_star_z;
                    end
                    break;
                end
            else
            % Step 6: Find the minimum of uncover value.
                min_Num = min(min(c_matrix_p(~cover_row,~cover_col)));
                % add min to cover row
                c_matrix_p(cover_row,:)=c_matrix_p(cover_row,:) + min_Num;
                % subtract min from uncover column
                c_matrix_p(:,~cover_col)=c_matrix_p(:,~cover_col) - min_Num;
            end
        end

    end
end