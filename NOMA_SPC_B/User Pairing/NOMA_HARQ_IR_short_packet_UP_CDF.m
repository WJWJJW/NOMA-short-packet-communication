clc; clear variables; close all;
N = 1e6; % number of channel tap
ncluster = 5:20;  % number of cluster (number of user  = 2K)
NN = 256; % number of information bit
N1 = NN;
N2 = NN;

NNN = 10000; % number of Monte Carlo

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

sum_En_UPG_opt_M_j = zeros(NNN,length(Pt));
sum_En_HAP_opt_M_j = zeros(NNN,length(Pt));

parfor u=1:length(ncluster)
    for jj=1:NNN
        h = (randn(1,N)+1i*randn(1,N));
        lamda = mean(abs(h).^2);

        K = ncluster(u);
        user_distance = randi([50 300],1,2*K);
%         target_BLER = (1e-4 - 1e-5).*rand(1,2*K) + 1e-5;
        exponent = 4 + 4*rand(1, 2*K);
        target_BLER = 10.^-exponent;
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
        
        % Enhanced Hungarian Algorithm Pairing
        [sum_En_HAP_opt_M_j(jj,u), ~, ~] =...
            En_HAP(user_distance, NN, K, target_BLER, rho, eta, lamda);
        
        % Enhanced User Pre-Grouping
        [sum_En_UPG_opt_M_j(jj,u), ~, ~] =...
            En_UPG_opt_delta(user_distance, NN, K, target_BLER, rho, eta, lamda);
        
    end
end


sum_RP_opt_M = mean(sum_RP_opt_M_j);
sum_UPG_opt_M = mean(sum_UPG_opt_M_j);
sum_HAP_opt_M = mean(sum_HAP_opt_M_j);
sum_OMA_opt_M = mean(sum_OMA_opt_M_j);
sum_NLUPA_opt_M = mean(sum_NLUPA_opt_M_j);
sum_SAP_opt_M = mean(sum_SAP_opt_M_j);

sum_En_UPG_opt_M = mean(sum_En_UPG_opt_M_j);
sum_En_HAP_opt_M = mean(sum_En_HAP_opt_M_j);

total_resource = 1.5e7;

RP_CDF = min(sum_RP_opt_M/total_resource,1);
UPG_CDF = min(sum_UPG_opt_M/total_resource,1);
HAP_CDF = min(sum_HAP_opt_M/total_resource,1);
OMA_CDF = min(sum_OMA_opt_M/total_resource,1);
NLUPA_CDF = min(sum_NLUPA_opt_M/total_resource,1);
SAP_CDF = min(sum_SAP_opt_M/total_resource,1);

En_UPG_CDF = min(sum_En_UPG_opt_M/total_resource,1);
En_HAP_CDF = min(sum_En_HAP_opt_M/total_resource,1);

figure (1)

plot(ncluster, sum_UPG_opt_M,'o','Color',[1 0.5 0]);
hold on; grid on;
plot(ncluster, sum_NLUPA_opt_M, 'm');
plot(ncluster, sum_HAP_opt_M, '.g');
plot(ncluster, sum_RP_opt_M,'b');
plot(ncluster, sum_OMA_opt_M,'c');
plot(ncluster, sum_SAP_opt_M,'r');

plot(ncluster, sum_En_UPG_opt_M, '--', 'Color',[1 0.5 0]);
plot(ncluster, sum_En_HAP_opt_M, '--g');

xlabel('Number of Cluster');
ylabel('Blocklength (Channel uses)');
legend('User Pre-Grouping', 'User Pre-Grouping NLUPA',...
    'Hungarian Algorithm Pairing', 'Random Pairing', 'OMA',...
    'Simulated Annealing Pairing', 'Enhanced User Pre-Grouping',...
    'Enhanced Hungarian Algorithm Pairing');

set(gca, 'FontName', 'Times New Roman');

figure (2)

plot(10:2:40, RP_CDF,'--b', 'linewidth', 1.5);
hold on; grid on;
plot(10:2:40, UPG_CDF, '-o', 'Color',[1 0.5 0], 'linewidth', 1.5);
plot(10:2:40, NLUPA_CDF, 'm', 'linewidth', 1.5);
plot(10:2:40, HAP_CDF, 'g', 'linewidth', 1.5);
plot(10:2:40, SAP_CDF, '-s', 'Color', [0.3010 0.7450 0.9330], 'linewidth', 1.5);
plot(10:2:40,OMA_CDF,'c', 'linewidth', 1.5);

plot(10:2:40, En_UPG_CDF, '--', 'Color',[1 0.5 0], 'linewidth', 1.5);
plot(10:2:40, En_HAP_CDF, '--g', 'linewidth', 1.5);


xlabel('Number of User');
ylabel('CDF');
legend('RP', 'UPG w/o RGP ', 'NLUPA',...
    'HAP', 'SAP', 'OMA',...
    'UPG','HRGP');

set(gca, 'FontName', 'Times New Roman');


% create a new pair of axes inside current figure
axes('position',[.65 .175 .25 .25])
box on % put box around new pair of axes
indexOfInterest = 30:2:34;

plot(indexOfInterest, SAP_CDF(11:13), '-s', 'Color', [0.3010 0.7450 0.9330], 'linewidth', 1.5);
hold on; grid on;
plot(indexOfInterest, En_HAP_CDF(11:13), '--g', 'linewidth', 1.5);


