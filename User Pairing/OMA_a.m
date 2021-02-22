% OMA
% assign reliability constraint fixedly
function [sum_opt_M, opt_M] = OMA_a(user_distance, N, K, eplsion1R, eplsion2R, rho, beta, OMA_PA, eta, lamda)
lamda_d = (1/2.*user_distance.^-eta).*lamda;
rho1 = OMA_PA/beta * rho;
rho2 = (1-OMA_PA)/(1-beta) * rho;


OMA_opt_M1 = N./(beta.*log2(1+(rho1.*eplsion1R.*lamda_d(1:K))));
OMA_opt_M2 = N./((1-beta).*log2(1+(rho2.*eplsion2R.*lamda_d(K+1:2*K))));

sum_opt_M = sum(OMA_opt_M1)+sum(OMA_opt_M2);
opt_M = [OMA_opt_M1, OMA_opt_M2];
end