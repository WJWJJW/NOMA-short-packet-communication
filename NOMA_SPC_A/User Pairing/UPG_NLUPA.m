% User Pre-Grouping (UPG) (From optimal UP paper NLUPA)
function [sum_opt_M, opt_M, User_pre_grouping] = UPG_NLUPA(user_distance, N, K, eplsion1R, eplsion2R, rho, eta, lamda)
    User_pre_grouping = zeros(K,2);
    % Paring processing
    for ii=1:K
        User_pre_grouping(ii,1) = user_distance(ii);
        User_pre_grouping(ii,2) = user_distance(2*K +1-ii);
    end
    
     % Total blocklength for User Pre-Grouping
    [sum_opt_M, opt_M] = M_cal_Mod(N,User_pre_grouping,K,eplsion1R,eplsion2R,rho,eta,lamda);
end