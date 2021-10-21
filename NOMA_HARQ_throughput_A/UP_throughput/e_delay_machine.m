% This function calcultate the expeted delay of user
function [expected_delay]=e_delay_machine(N1, T, user_distance, pair_idx, K,rho,eta,target_BLER, lamda)
epsilon2 = zeros(T, T+1);
epsilon1 = zeros(T, T+1);

epsilon2(:,1) = 1;
epsilon1(:,1) = 1;
opt_a1 = zeros(1,K);
opt_M = zeros(1,K);
T_cu_2 = zeros(K,T);
T_cu_1 = zeros(K,T);
expected_delay = zeros(K,2,T,'double');
user_pairing = user_distance(pair_idx);
target_BLER = target_BLER(pair_idx);
Delay = 200;
%Calculate optimal alpha1 and optimal blocklength
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
            opt_a1(ii) = (-lamda_d(1)*eplsion1R + sqrt((lamda_d(1)*eplsion1R)^2+...
                                4*(lamda_d(2)*eplsion2R_m)^2*rho*(lamda_d(1)*eplsion1R-lamda_d(2)*eplsion2R_m))) ...
                                    / (2*lamda_d(2)*eplsion2R_m*rho*(lamda_d(1)*eplsion1R-lamda_d(2)*eplsion2R_m));
            opt_M(ii) = N1/log2(1+((-lamda_d(1)*eplsion1R + sqrt((lamda_d(1)*eplsion1R)^2+4*(lamda_d(2)*eplsion2R_m)^2*rho*(lamda_d(1)*eplsion1R-lamda_d(2)*eplsion2R_m)))...
                                 / (2*lamda_d(2)*eplsion2R_m)));
                             
            for tt = 1:T
                sub_M = opt_M(ii) / tt;
                inc_M = opt_M(ii) / tt;
                for mm = 1:tt
                    % Average BLER (HARQ-IR)
                    w1 = 2^(N1/sub_M)-1;
                    w2 = 2^(N1/sub_M)-1;
                    epsilon2_high_SNR = w2/(lamda_d(2)*rho*(1-opt_a1(ii)-opt_a1(ii)*w2));
                    epsilon12_high_SNR = w2/(lamda_d(1)*rho*(1-opt_a1(ii)-opt_a1(ii)*w1));
                    epsilon11_high_SNR = w1/(lamda_d(1)*rho*opt_a1(ii));
                    
                    epsilon2(tt,mm+1) = epsilon2_high_SNR;
                    epsilon1(tt,mm+1) = epsilon12_high_SNR+epsilon11_high_SNR;

                    sub_M = sub_M + inc_M;
                end
                % Expected delay
                if tt == 1
                    T_cu_2(ii,tt) = opt_M(ii)*epsilon2(tt,1);
                    T_cu_1(ii,tt) = opt_M(ii)*epsilon1(tt,1);
                else
                    T_cu_2(ii,tt) = sum(inc_M*epsilon2(tt,1:tt))+...
                                            Delay.*sum(epsilon2(tt,1:tt));
                    T_cu_1(ii,tt) = sum(inc_M*epsilon1(tt,1:tt))+...
                                            Delay.*sum(epsilon1(tt,1:tt));  
                end
                expected_delay(ii,1,tt) = T_cu_1(ii,tt);
                expected_delay(ii,2,tt) = T_cu_2(ii,tt);
%                 N_expected_2(tt) = N2*(1-epsilon2(tt,mm+1));
%                 N_expected_1(tt) = N1*(1-epsilon1(tt,mm+1));
            end
%             % Calculate throughput
%             throughput_2 = N_expected_2./ T_cu_2;
%             throughput_1 = N_expected_1./ T_cu_1;
    end
    
end