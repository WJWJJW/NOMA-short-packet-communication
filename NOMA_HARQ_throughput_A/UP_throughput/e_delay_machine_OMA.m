% This function calcultate the expeted delay of user
function [expected_delay]=e_delay_machine_OMA(N1, T, user_distance, K, rho,eta, beta, OMA_PA, target_BLER, lamda)
epsilon = zeros(T, T+1);
epsilon(:,1) = 1;

T_cu = zeros(2*K,T);

expected_delay = zeros(2*K,T,'double');

Delay = 200;
%Calculate optimal alpha1 and optimal blocklength
rho1 = OMA_PA/beta * rho;
OMA_opt_M = zeros(1,2*K,'double');

for ii=1:2*K
    lamda_d = (1/2*user_distance(ii)^-eta)*lamda;
    OMA_opt_M(ii) = N1 / (beta*log2(1+(rho1*target_BLER(ii)*lamda_d)));
    for tt = 1:T
        sub_M = OMA_opt_M(ii) / tt;
        inc_M = OMA_opt_M(ii) / tt;
        for mm = 1:tt
            % Average BLER (HARQ-IR)
            w = 2^(N1/sub_M)-1;
            epsilon_high_SNR = w/(lamda_d*rho1);
            epsilon(tt,mm+1) = epsilon_high_SNR;
            sub_M = sub_M + inc_M;
        end
        % Expected delay
        if tt == 1
            T_cu(ii,tt) = OMA_opt_M(ii)*epsilon(tt,1);
        else
            T_cu(ii,tt) = sum(inc_M*epsilon(tt,1:tt))+...
                                    Delay.*sum(epsilon(tt,1:tt));
        end
        expected_delay(ii,tt) = T_cu(ii,tt);
    end
end

end