clc; clear variables; close all;
N = 1e6; % number of channel tap
ncluster = 5:20;  % number of cluster (number of user  = 2K)
NN = 80; % number of information bit
N1 = NN;
N2 = NN;

NNN = 10000; % number of Monte Carlo

eplsion1R = 10^-5;
eplsion2R = 10^-4;
delta = 1/2;

Pt = 20;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10^(No/10);   %Noise power (linear scale)

rho = pt/ no;
RHO = pow2db(rho);

eta = 4;

sum_RP_opt_M_j = zeros(NNN,length(ncluster));
sum_UPG_opt_M_j = zeros(NNN,length(ncluster));
sum_HAP_opt_M_j = zeros(NNN,length(ncluster));


sum_RP_opt_M = zeros(1,length(ncluster));
sum_UPG_opt_M = zeros(1,length(ncluster));
sum_HAP_opt_M = zeros(1,length(ncluster));



for u=1:length(ncluster)
    for jj=1:NNN
        h = (randn(1,N)+1i*randn(1,N));
        lamda = mean(abs(h).^2);

        K = ncluster(u);
        user_distance = randi([50 330],1,2*K);
%         if u == 1
%             user_distance = randi([10 330],1,2*K);
%         else
%             user_distance_add = randi([10 330],1,2);
%             user_distance = [user_distance user_distance_add];
%         end

        user_distance = sort(user_distance);

        % User Pre-Grouping
        [sum_UPG_opt_M_j(jj,u), ~, U] =...
            UPG_opt_delta(user_distance, NN, K, eplsion1R, eplsion2R, rho, eta, lamda, delta);
        
        
        % Random Paring (RP)
        [sum_RP_opt_M_j(jj,u), ~, ~]=...
            RP(user_distance, NN, K, eplsion1R, eplsion2R, rho, eta, lamda, delta);


        % Hungarian Algorithm Pairing
        [sum_HAP_opt_M_j(jj,u), ~, H] =...
            HAP(user_distance, NN, K, eplsion1R, eplsion2R, rho, eta, lamda, delta);
        
%         U
%         sum_UPG_opt_M_j(jj,u)
%         H
%         sum_HAP_opt_M_j(jj,u)
       
    end
end


sum_RP_opt_M = mean(sum_RP_opt_M_j);
sum_UPG_opt_M = mean(sum_UPG_opt_M_j);
sum_HAP_opt_M = mean(sum_HAP_opt_M_j);


figure (1)

plot(ncluster, sum_UPG_opt_M,'og');
hold on; grid on;
plot(ncluster, sum_HAP_opt_M,'Color',[1 0.5 0]);
plot(ncluster, sum_RP_opt_M,'r');


ylabel('blocklength');
legend('User Pre-Grouping', 'Hungarian Algorithm Pairing', 'Random Pairing');