function [sum_opt_M, cur_opt_M] = M_cal(N1,user_pairing,K,eplsion1R,eplsion2R,rho,eta, lamda, s)
cur_opt_M = zeros(1,K,'double');
    for ii =1:K       
        lamda_d = (1/2.*user_pairing(ii,:).^-eta).*lamda;
        
        contraint = (lamda_d(1)*eplsion1R)/(lamda_d(2)*eplsion2R);
        theta = contraint/s;
        
        check = (lamda_d(2)*eplsion2R*rho + 2)*lamda_d(1)*eplsion1R - (lamda_d(2)*eplsion2R.*rho + 4)*lamda_d(2)*eplsion2R > 0;
        switch (check)
            case 0
                eplsion2R_m = theta*eplsion2R;
            case 1
                eplsion2R_m = eplsion2R;
        end
        cur_opt_M(ii) = N1/log2(1+((-lamda_d(1)*eplsion1R + ...
                        sqrt((lamda_d(1)*eplsion1R)^2+4*(lamda_d(2)*eplsion2R_m)^2*rho*(lamda_d(1)*eplsion1R-lamda_d(2)*eplsion2R_m)))...
                        / (2*lamda_d(2)*eplsion2R_m)));          
    end
    sum_opt_M = sum(cur_opt_M(:));
end