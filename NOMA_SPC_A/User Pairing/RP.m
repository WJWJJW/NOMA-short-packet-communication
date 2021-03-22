% Random Paring (RP)
function [sum_opt_M, opt_M, RP_user_pairing]=RP(user_distance, N, K, target_BLER, rho, eta, lamda)
    RP_indices = randperm(2*K);
    RP_indices = sort(reshape(RP_indices,K,2),2);
    RP_user_pairing = user_distance(RP_indices);
    target_BLER_pair = target_BLER(RP_indices);

    
    % Total blocklength for random pairing
    [sum_opt_M, opt_M] = M_cal_Mod(N,RP_user_pairing, K, target_BLER_pair,rho,eta,lamda);
end