function [thred1_checker, thred2_checker] = thred_checker(user_pairing, K, eplsion1R, eplsion2R, rho, eta, lamda)
    dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);
    thred1_checker = zeros(K,1);
    thred2_checker = zeros(K,1);
    for ii=1:K
        lamda_d = (1/2.*user_pairing(ii,:).^-eta).*lamda;
        thred1_checker(ii) = double((user_pairing(ii,2) / user_pairing(ii,1)) > dis_thred1);

        check = (lamda_d(2)*eplsion2R*rho + 2)*lamda_d(1)*eplsion1R - (lamda_d(2)*eplsion2R.*rho + 4)*lamda_d(2)*eplsion2R > 0;
        switch (check)
            case 0
                thred2_checker(ii) = 0;
            case 1
                thred2_checker(ii) = 1;
        end


    end

end