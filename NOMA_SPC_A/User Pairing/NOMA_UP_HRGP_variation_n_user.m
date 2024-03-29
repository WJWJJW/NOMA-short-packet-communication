clc; clear variables; close all;
N = 1e6; % number of channel tap
ncluster = 5:40;  % number of cluster (number of user  = 2K)
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

sum_EP_opt_M_j = zeros(NNN,length(ncluster));
sum_RP_opt_M_j = zeros(NNN,length(ncluster));
sum_UPG_opt_M_j = zeros(NNN,length(ncluster));
sum_HAP_opt_M_j = zeros(NNN,length(ncluster));
sum_OMA_opt_M_j = zeros(NNN,length(ncluster));
sum_NLUPA_opt_M_j = zeros(NNN,length(ncluster));
sum_SAP_opt_M_j = zeros(NNN,length(ncluster));


sum_En_HAP_opt_M_j = zeros(NNN,length(ncluster));

sum_En_HRGP_max_d_opt_M_j = zeros(NNN,length(ncluster));
sum_En_HRGP_min_d_opt_M_j = zeros(NNN,length(ncluster));
sum_En_HRGP_min_BL_opt_M_j = zeros(NNN,length(ncluster));
sum_En_HRGP_min_BL_2_opt_M_j = zeros(NNN,length(ncluster));

parfor u=1:length(ncluster)
    K = ncluster(u);
    for jj=1:NNN
        h = (randn(1,N)+1i*randn(1,N));
        lamda = mean(abs(h).^2);

        user_distance = randi([50 300],1,2*K);
%         target_BLER = (1e-4 - 1e-5).*rand(1,2*K) + 1e-5;
        target_BLER = zeros(2*K, 1);
        target_BLER(1:K,:) = 1e-7;
        target_BLER(K+1:2*K,:) = 1e-4;
        
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
        
        % Enhanced HRGP max d
        [sum_En_HRGP_max_d_opt_M_j(jj,u), ~, ~] =...
            En_HRGP_max_d(user_distance, NN, K, target_BLER, rho, eta, lamda);
        
        % Enhanced HRGP min d
        [sum_En_HRGP_min_d_opt_M_j(jj,u), ~, ~] =...
            En_HRGP_min_d(user_distance, NN, K, target_BLER, rho, eta, lamda);
        
        % Enhanced HRGP min Blocklength
        [sum_En_HRGP_min_BL_opt_M_j(jj,u), ~, ~] =...
            En_HRGP_min_BL(user_distance, NN, K, target_BLER, rho, eta, lamda);
        
         % Enhanced HRGP min Blocklength
        [sum_En_HRGP_min_BL_2_opt_M_j(jj,u), ~, ~] =...
            En_HRGP_min_BL2(user_distance, NN, K, target_BLER, rho, eta, lamda);
               
    end
    
end

% sum_EP_opt_M = mean(sum_EP_opt_M_j(:,1:4));
sum_RP_opt_M = mean(sum_RP_opt_M_j);
sum_HAP_opt_M = mean(sum_HAP_opt_M_j);
sum_OMA_opt_M = mean(sum_OMA_opt_M_j);
sum_NLUPA_opt_M = mean(sum_NLUPA_opt_M_j);
sum_SAP_opt_M = mean(sum_SAP_opt_M_j);
sum_En_HAP_opt_M = mean(sum_En_HAP_opt_M_j);

sum_En_HRGP_max_d_opt_M = mean(sum_En_HRGP_max_d_opt_M_j);
sum_En_HRGP_min_d_opt_M = mean(sum_En_HRGP_min_d_opt_M_j);
sum_En_HRGP_min_BL_opt_M = mean(sum_En_HRGP_min_BL_opt_M_j);
sum_En_HRGP_min_BL_2_opt_M = mean(sum_En_HRGP_min_BL_2_opt_M_j);


figure (1)

plot(2*ncluster(1):2:2*ncluster(end), sum_NLUPA_opt_M, 'm', 'linewidth', 1.5);
hold on; grid on;
plot(2*ncluster(1):2:2*ncluster(end), sum_HAP_opt_M, '.g', 'linewidth', 1.5);
plot(2*ncluster(1):2:2*ncluster(end), sum_RP_opt_M,'b', 'linewidth', 1.5);
plot(2*ncluster(1):2:2*ncluster(end), sum_OMA_opt_M,'c', 'linewidth', 1.5);
plot(2*ncluster(1):2:2*ncluster(end), sum_SAP_opt_M,'-s', 'Color', [0.3010 0.7450 0.9330], 'linewidth', 1.5);
plot(2*ncluster(1):2:2*ncluster(end), sum_En_HAP_opt_M, '--g', 'linewidth', 1.5);

plot(2*ncluster(1):2:2*ncluster(end), sum_En_HRGP_max_d_opt_M, '-og', 'linewidth', 1.5);
plot(2*ncluster(1):2:2*ncluster(end), sum_En_HRGP_min_d_opt_M, '-sg', 'linewidth', 1.5);

plot(2*ncluster(1):2:2*ncluster(end), sum_En_HRGP_min_BL_opt_M, '-*g', 'linewidth', 1.5, 'MarkerEdgeColor','b');
plot(2*ncluster(1):2:2*ncluster(end), sum_En_HRGP_min_BL_2_opt_M, '-xg', 'linewidth', 1.5, 'MarkerEdgeColor','b');

xlabel('Number of Cluster');
ylabel('Blocklength (Channel uses)');
legend('User Pre-Grouping NLUPA',...
    'Hungarian Algorithm Pairing', 'Random Pairing', 'OMA',...
    'Simulated Annealing Pairing',...
    'Enhanced Hungarian Algorithm Pairing','HRGP max d','HRGP min d',...
        'HRGP min Blocklength','HRGP min Blocklength 2');

set(gca, 'FontName', 'Times New Roman');

figure (2)


plot(2*ncluster(1):2:2*ncluster(end), sum_En_HAP_opt_M, '--g', 'linewidth', 1.5);
hold on; grid on;
plot(2*ncluster(1):2:2*ncluster(end), sum_En_HRGP_max_d_opt_M, '-og', 'linewidth', 1.5);
plot(2*ncluster(1):2:2*ncluster(end), sum_En_HRGP_min_d_opt_M, '-sg', 'linewidth', 1.5);

plot(2*ncluster(1):2:2*ncluster(end), sum_En_HRGP_min_BL_opt_M, '-*g', 'linewidth', 1.5, 'MarkerEdgeColor','b');
plot(2*ncluster(1):2:2*ncluster(end), sum_En_HRGP_min_BL_2_opt_M, '-xg', 'linewidth', 1.5, 'MarkerEdgeColor','b');

xlabel('Number of Cluster');
ylabel('Blocklength (Channel uses)');
legend('HRGP','HRGP max d','HRGP min d',...
        'HRGP min Blocklength','HRGP min Blocklength 2');

set(gca, 'FontName', 'Times New Roman');


