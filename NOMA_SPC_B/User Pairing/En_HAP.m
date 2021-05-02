% Enhanced Hungarian Algorithm Pairing 
function [sum_opt_M, opt_M, Hungarian_pairing]=En_HAP(user_distance, N, K, target_BLER, rho, eta, lamda)
    user_idx = 1:2*K;
    UPG_idx = zeros(K,2);
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
    User_pre_grouping = user_distance(UPG_idx);
    target_BLER_pair = target_BLER(UPG_idx);

    c_matrix = zeros(K,K);
    
    % Generate cost matrix consist of blocklength for all possible pair
    for ii=1:K
        for jj=1:K
            if User_pre_grouping(ii,1) < User_pre_grouping(jj,2)
                [c_matrix(ii,jj)] = M_cal_Mod(N,[User_pre_grouping(ii,1),User_pre_grouping(jj,2)],1,...
                            [target_BLER_pair(ii,1),target_BLER_pair(jj,2)],rho,eta,lamda);
            else 
                [c_matrix(ii,jj)] = M_cal_Mod(N,[User_pre_grouping(jj,2),User_pre_grouping(ii,1)],1,...
                            [target_BLER_pair(jj,2),target_BLER_pair(ii,1)],rho,eta,lamda);
            end
        end
    end
    
    % Apply Hungarian Algorithm
    
    [starZ] = Hungarian_algorithm(c_matrix,K);
    [near,far] = find(starZ);
    Hungarian_pairing(:,1) = User_pre_grouping(near,1);
    Hungarian_pairing(:,2) = User_pre_grouping(far,2);
    Hungarian_pairing = sort(Hungarian_pairing,2);
    sum_opt_M = sum(c_matrix(starZ == 1));
    opt_M = c_matrix(starZ == 1)';

end
