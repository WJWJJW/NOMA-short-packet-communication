function [thred1_checker, thred2_checker] = thred_checker(user_pairing, K, eplsion1R, eplsion2R, rho, eta)
    dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);
    thred1_checker = zeros(K,1);
    thred2_checker = zeros(K,1);
    for ii=1:K
        thred1_checker(ii) = double((user_pairing(ii,2) / user_pairing(ii,1)) > dis_thred1);
        if (user_pairing(ii,2)^-eta*eplsion2R*rho + 2)*user_pairing(ii,1)^-eta*eplsion1R... 
            - (user_pairing(ii,2)^-eta*eplsion2R.*rho + 4)*user_pairing(ii,2)^-eta*eplsion2R > 0
            thred2_checker(ii) = 1;
        end
    end

end