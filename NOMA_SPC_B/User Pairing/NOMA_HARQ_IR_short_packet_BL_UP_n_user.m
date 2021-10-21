clc; clear variables; close all;
N = 1e6; % number of channel tap
ncluster = 5:20;  % number of cluster (number of user  = 2K)
NN = 256; % number of information bit
N1 = NN;
N2 = NN;

NNN = 10000; % number of Monte Carlo

eplsion1R = 10^-5;
eplsion2R = 10^-4;


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
    K = ncluster(u);
    for jj=1:NNN
        h = (randn(1,N)+1i*randn(1,N));
        lamda = mean(abs(h).^2);

        user_distance = randi([50 300],1,2*K);
        target_BLER = (1e-4 - 1e-5).*rand(1,2*K) + 1e-5;

        
        user_distance = sort(user_distance);
        
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
        
        % OMA
        [sum_OMA_opt_M_j(jj,u),~] = ...
            OMA(user_distance, NN, K, target_BLER, rho, beta, OMA_PA, eta, lamda);
        
        % Enhanced Hungarian Algorithm Pairing
        [sum_En_HAP_opt_M_j(jj,u), ~, ~] =...
            En_HAP(user_distance, NN, K, target_BLER, rho, eta, lamda);
        
    end
end


sum_RP_opt_M = mean(sum_RP_opt_M_j);
sum_HAP_opt_M = mean(sum_HAP_opt_M_j);
sum_OMA_opt_M = mean(sum_OMA_opt_M_j);
sum_NLUPA_opt_M = mean(sum_NLUPA_opt_M_j);
sum_SAP_opt_M = mean(sum_SAP_opt_M_j);
sum_En_HAP_opt_M = mean(sum_En_HAP_opt_M_j);

figure (1)

plot(ncluster, sum_NLUPA_opt_M, 'm');
hold on; grid on;
plot(ncluster, sum_HAP_opt_M, '.g');
plot(ncluster, sum_RP_opt_M,'b');
plot(ncluster, sum_OMA_opt_M,'c');
plot(ncluster, sum_SAP_opt_M,'-s', 'Color', [0.3010 0.7450 0.9330]);
plot(ncluster, sum_En_HAP_opt_M, '--g');


xlabel('Number of Cluster');
ylabel('Blocklength (Channel uses)');
legend('User Pre-Grouping NLUPA',...
    'Hungarian Algorithm Pairing', 'Random Pairing', 'OMA',...
    'Simulated Annealing Pairing',...
    'Enhanced Hungarian Algorithm Pairing');

set(gca, 'FontName', 'Times New Roman');

figure (2)
plot(ncluster, sum_HAP_opt_M, '.g');
hold on; grid on;
plot(ncluster, sum_SAP_opt_M,'-s', 'Color', [0.3010 0.7450 0.9330]);
plot(ncluster, sum_En_HAP_opt_M, '--g');


xlabel('Number of Cluster');
ylabel('Blocklength (Channel uses)');
legend('Hungarian Algorithm Pairing',...
    'Simulated Annealing Pairing',...
    'Enhanced Hungarian Algorithm Pairing');

set(gca, 'FontName', 'Times New Roman');



