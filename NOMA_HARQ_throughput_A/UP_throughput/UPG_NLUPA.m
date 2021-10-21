% User Pre-Grouping (UPG) (From optimal UP paper NLUPA)
function [sum_opt_M, opt_M, User_pre_grouping, UPG_idx] = UPG_NLUPA(user_distance, N, K, target_BLER, rho, eta, lamda)
    user_idx = 1:2*K;
    UPG_idx = zeros(K,2);

    % Paring processing
    for ii=1:K
        UPG_idx(ii,1) = user_idx(ii);
        UPG_idx(ii,2) = user_idx(2*K +1-ii);
    end
    User_pre_grouping = user_distance(UPG_idx);
    target_BLER_pair = target_BLER(UPG_idx);
    
     % Total blocklength for User Pre-Grouping
    [sum_opt_M, opt_M] = M_cal_Mod(N,User_pre_grouping,K,target_BLER_pair,rho,eta,lamda);
end