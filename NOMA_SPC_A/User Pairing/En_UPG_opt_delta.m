% Enhanced User Pre-Grouping optimal delta (UPG_opt_delta)
function [sum_opt_M, opt_M, User_pre_grouping] = En_UPG_opt_delta(user_distance, N, K, target_BLER, rho, eta, lamda)

    user_idx = 1:2*K;
    UPG_idx = zeros(K,2);
%     NU = zeros(K,1);
%     FU = zeros(K,1);
    re_pairing_idx_check = zeros(K,1);
    % Paring processing
    for ii = 1:K
        UPG_idx(ii,1) = user_idx(ii);
        UPG_idx(ii,2) = user_idx(K + ii);
        
        d_min = min_d_finder(user_distance(ii), target_BLER(ii), target_BLER(K + ii), rho);
        if d_min > user_distance(K + ii)
            re_pairing_idx_check(ii) = 1;
        end
    end
    % Regrouping Process
    if size(re_pairing_idx_check(re_pairing_idx_check == 1),1) == 1
        re_pairing_idx_check(1) = 1;
    end
    re_pair_user_idx = user_idx(UPG_idx(re_pairing_idx_check==1,:));
    t = 1:2:numel(re_pair_user_idx);
    for ii = 1: size(re_pair_user_idx,1)
        UPG_idx(re_pair_user_idx(ii),1) = re_pair_user_idx(t(ii));
        UPG_idx(re_pair_user_idx(ii),2) = re_pair_user_idx(t(ii)+1);
    end
%     NU = sort(UPG_idx(:,1));
%     FU = sort(UPG_idx(:,2));
%     for ii=1:K
%         UPG_idx(ii,1) = NU(ii);
%         UPG_idx(ii,2) = FU(ii);
%     end
    User_pre_grouping = user_distance(UPG_idx);
    target_BLER_pair = target_BLER(UPG_idx);
    
     % Total blocklength for User Pre-Grouping
    [sum_opt_M, opt_M] = M_cal_Mod(N,User_pre_grouping,K,target_BLER_pair,rho,eta,lamda);
end