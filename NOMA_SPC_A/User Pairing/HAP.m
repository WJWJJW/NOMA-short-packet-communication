% function [sum_opt_M, opt_M, Hungarian_pairing]=HAP(user_distance, N, K, eplsion1R, eplsion2R, rho, eta, lamda, delta)
%     n_set = user_distance(1:K);
%     f_set = user_distance(K+1:2*K);
%     c_matrix = zeros(K,K);
%     
%     for ii=1:K
%         for jj=1:K
% %             [c_matrix(ii,jj)] = M_cal(N,[n_set(ii),f_set(jj)],1,eplsion1R,eplsion2R,rho,eta,lamda,delta);
%             [c_matrix(ii,jj)] = M_cal_Mod(N,[n_set(ii),f_set(jj)],1,eplsion1R,eplsion2R,rho,eta,lamda);
%         end
%     end
%     ind_matrix = matchpairs(c_matrix, max(max(c_matrix)));
%     Hungarian_pairing(:,1) = n_set(ind_matrix(:,1));
%     Hungarian_pairing(:,2) = f_set(ind_matrix(:,2));
%     
%     sz = [K,K];
%     sum_opt_M = sum(c_matrix(sub2ind(sz,ind_matrix(:,1),ind_matrix(:,2))));
%     opt_M = c_matrix(sub2ind(sz,ind_matrix(:,1),ind_matrix(:,2)))';
% end
% 
function [sum_opt_M, opt_M, Hungarian_pairing]=HAP(user_distance, N, K, eplsion1R, eplsion2R, rho, eta, lamda)
    n_set = user_distance(1:K);
    f_set = user_distance(K+1:2*K);
    c_matrix = zeros(K,K);
    
    % Generate cost matrix consist of blocklength for all possible pair
    for ii=1:K
        for jj=1:K
%             [c_matrix(ii,jj)] = M_cal(N,[n_set(ii),f_set(jj)],1,eplsion1R,eplsion2R,rho,eta,lamda,delta);
            [c_matrix(ii,jj)] = M_cal_Mod(N,[n_set(ii),f_set(jj)],1,eplsion1R,eplsion2R,rho,eta,lamda);
        end
    end
    [starZ] = Hungarian_algorithm(c_matrix,K);
    [near,far] = find(starZ);
    Hungarian_pairing(:,1) = n_set(near);
    Hungarian_pairing(:,2) = f_set(far);
    sum_opt_M = sum(c_matrix(starZ == 1));
    opt_M = c_matrix(starZ == 1)';

end




% function [sum_opt_M, opt_M, Hungarian_pairing]=HAP(user_distance, N, K, eplsion1R, eplsion2R, rho, eta, lamda, delta)
%     n_set = user_distance(1:K);
%     f_set = user_distance(K+1:2*K);
%     c_matrix = zeros(K,K);
%     
%     for ii=1:K
%         for jj=1:K
% %             [c_matrix(ii,jj)] = M_cal(N,[n_set(ii),f_set(jj)],1,eplsion1R,eplsion2R,rho,eta,lamda,delta);
%             [c_matrix(ii,jj)] = M_cal_Mod(N,[n_set(ii),f_set(jj)],1,eplsion1R,eplsion2R,rho,eta,lamda);
%         end
%     end
%     [sol_row, sol_col] = linear_sum_assignment(c_matrix);
%     Hungarian_pairing(:,1) = n_set(sol_row);
%     Hungarian_pairing(:,2) = f_set(sol_col);
%     
%     [sum_opt_M, opt_M] = M_cal_Mod(N,Hungarian_pairing,K,eplsion1R,eplsion2R,rho,eta,lamda);
% end














