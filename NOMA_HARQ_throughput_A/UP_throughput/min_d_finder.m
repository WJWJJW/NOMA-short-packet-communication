% Thm 4
function [min_d]=min_d_finder(d1,eplsion1R, eplsion2R, rho)
    eta = 4;   
    lamda1 = d1^-eta;
    beta1 = 0.5;
    beta2 = 0.5;
    rho1 = rho;
    rho2 = rho;
    [delta]=delta_finder(lamda1);
    term1 = beta1*log(1+eplsion1R*lamda1*rho1);
    term2 = log(1+((-1+sqrt(1+4*delta^2*rho*lamda1*eplsion1R*(1-delta))) / (2*delta)));
    min_d = ((exp((term1*term2)/((term1-term2)*beta2))-1) / (eplsion2R*rho2))^(-1/eta);
end