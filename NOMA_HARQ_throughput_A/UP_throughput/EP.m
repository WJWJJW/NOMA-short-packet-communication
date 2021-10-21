% Exhaustive Paring (EP)
function [sum_opt_M, opt_M, Exhaustive_pairing, EP_idx]=EP(user_distance, pair_idx, N, K, target_BLER, rho, eta, lamda)
    
    exhaustive_pairing  = user_distance(pair_idx);
    target_BLER_pair = target_BLER(pair_idx);
    [tmp_sum_M, tmp_M] = M_cal_Mod(N, exhaustive_pairing(:,:,1),K,target_BLER_pair(:,:,1),rho,eta,lamda);
    Exhaustive_pairing = exhaustive_pairing(:,:,1);
    EP_idx = pair_idx(:,:,1);
    
    for jj=2:length(exhaustive_pairing)
        % Total blocklength for exhaustive paring
        [sum_EP_M, EP_M] = M_cal_Mod(N, exhaustive_pairing(:,:,jj),K,target_BLER_pair(:,:,jj),rho,eta,lamda);
        
        if sum_EP_M < tmp_sum_M
            tmp_sum_M = sum_EP_M;
            tmp_M = EP_M;
            Exhaustive_pairing = exhaustive_pairing(:,:,jj);
            EP_idx = pair_idx(:,:,jj);
        end       
    end
    sum_opt_M = tmp_sum_M;
    opt_M = tmp_M;

end
