clc; clear variables; close all;
N = 1e6; % number of channel tap
ncluster = 5:20;  % number of cluster (number of user  = 2K)
NN = 256; % number of information bit
N1 = NN;
N2 = NN;

NNN = 1000; % number of Monte Carlo

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


sum_OMA_opt_M_j = zeros(NNN,length(ncluster));

sum_En_UPG_opt_M_j = zeros(NNN,length(Pt));
sum_En_HAP_opt_M_j = zeros(NNN,length(Pt));



parfor u=1:length(ncluster)
    for jj=1:NNN
        h = (randn(1,N)+1i*randn(1,N));
        lamda = mean(abs(h).^2);

        K = ncluster(u);
        user_distance = randi([50 300],1,2*K);
        
        exponent = 7 + zeros(1,K);
        target_BLER_f = 10.^-exponent;
        exponent = 4 + zeros(1,K);
        target_BLER_n = 10.^-exponent;

        target_BLER = [target_BLER_n target_BLER_f];
        
%         exponent = 6 + 2*rand(1,K);
%         target_BLER_n = 10.^-exponent;
%         exponent = 3 + 2*rand(1,K);
%         target_BLER_f = 10.^-exponent;
%         
%         target_BLER = [target_BLER_n target_BLER_f];
        
        user_distance = sort(user_distance);
      
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

sum_OMA_opt_M = mean(sum_OMA_opt_M_j);

sum_En_UPG_opt_M = mean(sum_En_UPG_opt_M_j);
sum_En_HAP_opt_M = mean(sum_En_HAP_opt_M_j);


% figure (1)
% bar(ncluster, sum_En_HAP_opt_M, 'g');
% hold on; grid on;
% bar(ncluster, sum_En_UPG_opt_M,'r');
% 
% bar(ncluster, sum_OMA_opt_M,'b');
% 
% xlabel('Number of Cluster');
% ylabel('Blocklength (Channel uses)');
% legend('HAP','UPG', 'OMA');
% 
% set(gca, 'FontName', 'Times New Roman');

ret = [sum_En_HAP_opt_M;sum_En_UPG_opt_M;sum_OMA_opt_M];
figure (1)
bar(ncluster, ret');


xlabel('Number of Cluster');
ylabel('Blocklength (Channel uses)');
legend('HAP','UPG', 'OMA');

set(gca, 'FontName', 'Times New Roman');

