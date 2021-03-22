 function [sum_opt_M, opt_M, Hungarian_pairing]=HAP(user_distance, N, K, target_BLER, rho, eta, lamda)
    n_set = user_distance(1:K);
    f_set = user_distance(K+1:2*K);
    target_BLER_n = target_BLER(1:K);
    target_BLER_f = target_BLER(K+1:2*K);
    c_matrix = zeros(K,K);
    
    % Generate cost matrix consist of blocklength for all possible pair
    for ii=1:K
        for jj=1:K
%             [c_matrix(ii,jj)] = M_cal(N,[n_set(ii),f_set(jj)],1,eplsion1R,eplsion2R,rho,eta,lamda,delta);
            [c_matrix(ii,jj)] = M_cal_Mod(N,[n_set(ii),f_set(jj)],1,...
                            [target_BLER_n(ii),target_BLER_f(jj)],rho,eta,lamda);
        end
    end
    [starZ] = Hungarian_algorithm(c_matrix,K);
    [near,far] = find(starZ);
    Hungarian_pairing(:,1) = n_set(near);
    Hungarian_pairing(:,2) = f_set(far);
    sum_opt_M = sum(c_matrix(starZ == 1));
    opt_M = c_matrix(starZ == 1)';

end














