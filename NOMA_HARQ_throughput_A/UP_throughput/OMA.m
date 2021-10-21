% OMA
% assign reliability constraint according to pairing rule
function [sum_opt_M, OMA_opt_M] = OMA(user_distance, N, K, target_BLER, rho, beta, OMA_PA, eta, lamda)

rho1 = OMA_PA/beta * rho;
OMA_opt_M = zeros(1,K,'double');

for ii=1:2*K
    lamda_d = (1/2*user_distance(ii)^-eta)*lamda;
    OMA_opt_M(ii) = N / (beta*log2(1+(rho1*target_BLER(ii)*lamda_d)));
end

sum_opt_M = sum(OMA_opt_M);
end