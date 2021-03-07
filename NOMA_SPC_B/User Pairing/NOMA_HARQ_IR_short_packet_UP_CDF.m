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

Pt = 30;                    %Transmit Power in dBm
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
sum_NLUPA_opt_M_j = zeros(NNN,length(ncluster));
sum_SAP_opt_M_j = zeros(NNN,length(ncluster));

parfor u=1:length(ncluster)
    for jj=1:NNN
        h = (randn(1,N)+1i*randn(1,N));
        lamda = mean(abs(h).^2);

        K = ncluster(u);
        user_distance = randi([50 300],1,2*K);
        target_BLER = (1e-4 - 1e-5).*rand(1,2*K) + 1e-5;
%         if u == 1
%             user_distance = randi([10 330],1,2*K);
%         else
%             user_distance_add = randi([10 330],1,2);
%             user_distance = [user_distance user_distance_add];
%         end

        user_distance = sort(user_distance);

        % User Pre-Grouping
        [sum_UPG_opt_M_j(jj,u), ~, ~] =...
            UPG_opt_delta(user_distance, NN, K, target_BLER, rho, eta, lamda);
        
        % User Pre-Grouping NLUPA
        [sum_NLUPA_opt_M_j(jj,u), ~, ~] =...
            UPG_NLUPA(user_distance, NN, K, target_BLER, rho, eta, lamda);
        
        % Random Paring (RP)
        [sum_RP_opt_M_j(jj,u), ~, ~]=...
            RP(user_distance, NN, K, target_BLER, rho, eta, lamda);

        % Hungarian Algorithm Pairing
        [sum_HAP_opt_M_j(jj,u), ~, ~] =...
            HAP(user_distance, NN, K, target_BLER, rho, eta, lamda);
       
        % Simulated Annealing Pairing (SAP)
        [sum_SAP_opt_M_j(jj,u), ~, ~] =...
            SAP(user_distance, NN, K, target_BLER, rho, eta, lamda);
        
        % OMA (reliability constraint according to pair rule)
        [sum_OMA_opt_M_j(jj,u),~] = ...
            OMA(user_distance, NN, K, target_BLER, rho, beta, OMA_PA, eta, lamda);
        
    end
end


sum_RP_opt_M = mean(sum_RP_opt_M_j);
sum_UPG_opt_M = mean(sum_UPG_opt_M_j);
sum_HAP_opt_M = mean(sum_HAP_opt_M_j);
sum_OMA_opt_M = mean(sum_OMA_opt_M_j);
sum_NLUPA_opt_M = mean(sum_NLUPA_opt_M_j);
sum_SAP_opt_M = mean(sum_SAP_opt_M_j);


RP_CDF = min(sum_RP_opt_M/sum_HAP_opt_M(14),1);
UPG_CDF = min(sum_UPG_opt_M/sum_HAP_opt_M(14),1);
HAP_CDF = min(sum_HAP_opt_M/sum_HAP_opt_M(14),1);
OMA_CDF = min(sum_OMA_opt_M/sum_HAP_opt_M(14),1);
NLUPA_CDF = min(sum_NLUPA_opt_M/sum_HAP_opt_M(14),1);
SAP_CDF = min(sum_SAP_opt_M/sum_HAP_opt_M(14),1);

figure (1)

plot(ncluster, sum_UPG_opt_M,'Color',[1 0.5 0]);
hold on; grid on;
plot(ncluster, sum_NLUPA_opt_M, 'm');
plot(ncluster, sum_HAP_opt_M, 'g*');
plot(ncluster, sum_RP_opt_M,'b');
plot(ncluster, sum_OMA_opt_M,'c');
plot(ncluster, sum_SAP_opt_M,'r');

xlabel('Number of Cluster');
ylabel('Blocklength (Channel uses)');
legend('User Pre-Grouping', 'User Pre-Grouping NLUPA',...
    'Hungarian Algorithm Pairing', 'Random Pairing', 'OMA',...
    'Simulated Annealing Pairing');

set(gca, 'FontName', 'Times New Roman');

figure (2)

plot(10:2:40, RP_CDF,'--b');
hold on; grid on;
plot(10:2:40, UPG_CDF, '-o', 'Color',[1 0.5 0]);
plot(10:2:40, NLUPA_CDF, '-+m');
plot(10:2:40, HAP_CDF, '-.g');
plot(10:2:40, SAP_CDF, '-s', 'Color', [0.3010 0.7450 0.9330]);
plot(10:2:40,OMA_CDF,'c');


xlabel('Number of User');
ylabel('CDF');
legend('Random Pairing', 'User Pre-Grouping', 'User Pre-Grouping NLUPA',...
    'Hungarian Algorithm Pairing', 'Simulated Annealing Pairing', 'OMA');

set(gca, 'FontName', 'Times New Roman');
