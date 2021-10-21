% Enhanced HRGP max d
function [sum_opt_M, opt_M, Hungarian_pairing]=En_HRGP_max_d(user_distance, N, K, target_BLER, rho, eta, lamda)
    user_idx = 1:2*K;
    UPG_idx = zeros(K,2);
    re_pairing_idx_check = zeros(K,1);
    d_min = zeros(K,1);
    % Paring processing
    for ii = 1:K
        UPG_idx(ii,1) = user_idx(ii);
        UPG_idx(ii,2) = user_idx(K + ii);
        
        d_min(ii) = min_d_finder(user_distance(ii), target_BLER(ii), target_BLER(K + ii), rho);
        if d_min(ii) > user_distance(K + ii)
            re_pairing_idx_check(ii) = 1;
        end
    end
    % Regrouping Process
    % only 1 pair do not satisfy Thm. 4
    if size(re_pairing_idx_check(re_pairing_idx_check == 1),1) == 1
        re_pairing_idx_check(1) = 1;
    end
   
    
    if rem(numel(re_pairing_idx_check(re_pairing_idx_check==1)), 2) == 1 % |RG|/2 is odd
        [~, max_d_min_idx] = max(re_pairing_idx_check .* d_min);
        re_pairing_idx_check(max_d_min_idx) = 0;
    end % |RG|/2 is even, do nothing
    
    re_pair_user_idx = user_idx(UPG_idx(re_pairing_idx_check==1,:));
    
    NU = [reshape(UPG_idx(re_pairing_idx_check==0,1),1,numel(UPG_idx(re_pairing_idx_check==0,1))) ...
          reshape(re_pair_user_idx(1:2:size(re_pair_user_idx),:),1,numel(re_pair_user_idx(1:2:size(re_pair_user_idx),:)))];
    FU = [reshape(UPG_idx(re_pairing_idx_check==0,2),1,numel(UPG_idx(re_pairing_idx_check==0,2))) ...
          reshape(re_pair_user_idx(2:2:size(re_pair_user_idx),:),1,numel(re_pair_user_idx(2:2:size(re_pair_user_idx),:)))];  
%     NU = [reshape(UPG_idx(re_pairing_idx_check==0,1),1,numel(UPG_idx(re_pairing_idx_check==0,1))) ...
%           reshape(re_pair_user_idx(1:size(re_pair_user_idx)/2,:),1,numel(re_pair_user_idx(1:size(re_pair_user_idx)/2,:)))];
%     FU = [reshape(UPG_idx(re_pairing_idx_check==0,2),1,numel(UPG_idx(re_pairing_idx_check==0,2))) ...
%           reshape(re_pair_user_idx(size(re_pair_user_idx)/2+1:size(re_pair_user_idx),:),1,numel(re_pair_user_idx(size(re_pair_user_idx)/2+1:size(re_pair_user_idx),:)))];
    UPG_idx = [NU', FU'];
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
