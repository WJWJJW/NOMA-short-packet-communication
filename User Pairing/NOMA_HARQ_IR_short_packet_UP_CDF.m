clc; clear variables; close all;
N = 1e6; % number of channel tap
ncluster = 5:20;  % number of cluster (number of user  = 2K)
NN = 256; % number of information bit
N1 = NN;
N2 = NN;

NNN = 1000; % number of Monte Carlo

eplsion1R = 10^-5;
eplsion2R = 10^-4;
delta = 1/2;

Pt = 20;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

% AWGN
% BW = 10^7;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt/ no;
RHO = pow2db(rho);

beta = 0.5;
OMA_PA = 0.5;

eta = 4;

sum_RP_opt_M_j = zeros(NNN,length(ncluster));
sum_UPG_opt_M_j = zeros(NNN,length(ncluster));
sum_HAP_opt_M_j = zeros(NNN,length(ncluster));
sum_OMA_opt_M_j = zeros(NNN,length(ncluster));


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
        [sum_UPG_opt_M_j(jj,u), ~, UPG] =...
            UPG_opt_delta(user_distance, NN, K, eplsion1R, eplsion2R, rho, eta, lamda);
        
        
        % Random Paring (RP)
        [sum_RP_opt_M_j(jj,u), ~, ~]=...
            RP(user_distance, NN, K, eplsion1R, eplsion2R, rho, eta, lamda);

        % Hungarian Algorithm Pairing
        [sum_HAP_opt_M_j(jj,u), ~, ~] =...
            HAP(user_distance, NN, K, eplsion1R, eplsion2R, rho, eta, lamda);
       
        % OMA (reliability constraint according to pair rule)
        [sum_OMA_opt_M_j(jj,u),~] = ...
            OMA_b(UPG, NN, K, eplsion1R, eplsion2R, rho, beta, OMA_PA, eta, lamda);
        
    end
end


sum_RP_opt_M = mean(sum_RP_opt_M_j);
sum_UPG_opt_M = mean(sum_UPG_opt_M_j);
sum_HAP_opt_M = mean(sum_HAP_opt_M_j);
sum_OMA_opt_M = mean(sum_OMA_opt_M_j);


figure (1)

plot(ncluster, sum_UPG_opt_M,'og');
hold on; grid on;
plot(ncluster, sum_HAP_opt_M,'Color',[1 0.5 0]);
plot(ncluster, sum_RP_opt_M,'r');
plot(ncluster, sum_OMA_opt_M,'c');

xlabel('Number of cluster');
ylabel('Blocklength (Channel uses)');
legend('User Pre-Grouping', 'Hungarian Algorithm Pairing', 'Random Pairing', 'OMA');

set(gca, 'FontName', 'Times New Roman');