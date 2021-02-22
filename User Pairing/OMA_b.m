% OMA
% assign reliability constraint according to pairing rule
function [sum_opt_M, opt_M] = OMA_b(user_pairing, N, K, eplsion1R, eplsion2R, rho, beta, OMA_PA, eta, lamda)

rho1 = OMA_PA/beta * rho;
rho2 = (1-OMA_PA)/(1-beta) * rho;

OMA_opt_M1 = zeros(1,K,'double');
OMA_opt_M2 = zeros(1,K,'double');

for ii=1:K
    lamda_d = (1/2.*user_pairing(ii,:).^-eta).*lamda;
    OMA_opt_M1(ii) = N./(beta.*log2(1+(rho1.*eplsion1R.*lamda_d(1))));
    OMA_opt_M2(ii) = N./((1-beta).*log2(1+(rho2.*eplsion2R.*lamda_d(2))));
end

sum_opt_M = sum(OMA_opt_M1)+sum(OMA_opt_M2);
opt_M = [OMA_opt_M1, OMA_opt_M2];
end