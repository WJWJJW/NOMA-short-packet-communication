% This function calculate the blocklength by optimal delta via delta_finder
function [sum_opt_M, cur_opt_M] = M_cal_Mod(N1,user_pairing,K,target_BLER,rho,eta, lamda)
cur_opt_M = zeros(1,K,'double');
    for ii =1:K       
        lamda_d = (1/2.*user_pairing(ii,:).^-eta).*lamda;
        
        eplsion1R = target_BLER(ii,1);
        eplsion2R = target_BLER(ii,2);
        
        contraint = (lamda_d(1)*eplsion1R)/(lamda_d(2)*eplsion2R);
        [delta] = delta_finder (lamda_d(1)*rho*eplsion1R);
        theta = contraint*delta;
        
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