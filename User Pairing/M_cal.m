function [sum_opt_M, cur_opt_M] = M_cal(N1,user_pairing,K,eplsion1R,eplsion2R,rho,N,eta)
cur_opt_M = zeros(1,5);
    for ii =1:K
        h1 = sqrt(1/2*user_pairing(ii,1)^-eta)*(randn(1,N)+1i*randn(1,N));
        h2 = sqrt(1/2*user_pairing(ii,2)^-eta)*(randn(1,N)+1i*randn(1,N));
        
        lamda1 = mean(abs(h1).^2);
        lamda2 = mean(abs(h2).^2);
        
        contraint = (lamda1*eplsion1R)/(lamda2*eplsion2R);
        k = contraint/3;
        
        if user_pairing(ii,2)/user_pairing(ii,1) > (eplsion2R/eplsion1R)^0.25
%             oopt_a1(ii) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
%                       4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R))) ...
%                       / (2*lamda2*eplsion2R*rho*(lamda1*eplsion1R-lamda2*eplsion2R));
            cur_opt_M(ii) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R)))...
                                 / (2*lamda2*eplsion2R)));
        
        elseif user_pairing(ii,2)/user_pairing(ii,1) > (eplsion2R/eplsion1R)^0.25 && ...
               (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - (lamda2*eplsion2R.*rho + 4)*lamda2*eplsion2R < 0
        
%            oopt_a1(ii) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
%                   4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R))) ...
%                   / (2*k*lamda2*eplsion2R*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R));
           cur_opt_M(ii) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R)))...
                     / (2*k*lamda2*eplsion2R)));
    
        else
%             oopt_a1(ii) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
%                   4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R))) ...
%                   / (2*k*lamda2*eplsion2R*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R));
            cur_opt_M(ii) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R)))...
                     / (2*k*lamda2*eplsion2R)));
        end
                 
                 
    end
    sum_opt_M = sum(cur_opt_M(:));
end