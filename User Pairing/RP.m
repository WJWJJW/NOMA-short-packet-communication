% Random Paring (RP)
function [sum_opt_M, opt_M, RP_user_pairing]=RP(user_distance, N, K, eplsion1R, eplsion2R, rho, eta, lamda)
    RP_indices = randperm(2*K);
    RP_user_pairing = zeros(K,2);
    for ii=1:K
        RP_user_pairing(ii,:) = sort(user_distance(RP_indices(2*ii-1:2*ii)));
    end
    
    % Total blocklength for random pairing
    [sum_opt_M, opt_M] = M_cal_Mod(N,RP_user_pairing, K, eplsion1R,eplsion2R,rho,eta,lamda);
end