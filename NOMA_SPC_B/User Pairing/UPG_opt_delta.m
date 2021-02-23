% User Pre-Grouping optimal delta (UPG_opt_delta)
function [sum_opt_M, opt_M, User_pre_grouping] = UPG_opt_delta(user_distance, N, K, target_BLER, rho, eta, lamda)
    User_pre_grouping = zeros(K,2);
    target_BLER_pair = zeros(K,2);
    % Paring processing
    for ii=1:K
        User_pre_grouping(ii,1) = user_distance(ii);
        User_pre_grouping(ii,2) = user_distance(K + ii);
        target_BLER_pair(ii,1) = target_BLER(ii);
        target_BLER_pair(ii,2) = target_BLER(K + ii);
    end
    
     % Total blocklength for User Pre-Grouping
    [sum_opt_M, opt_M] = M_cal_Mod(N,User_pre_grouping,K,target_BLER_pair,rho,eta,lamda);
end