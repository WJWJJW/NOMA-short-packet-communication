clc; clear variables; close all;

N = 1e6; % number of channel tap
NNN = 1000; % number of Monte Carlo
K = 5  % number of cluster (number of user  = 2K)
NN = 256; % number of information bit
N1 = NN;
N2 = NN;


Pt = 20:2:30;               %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

% AWGN
% BW = 10^7;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt./ no;
RHO = pow2db(rho);

beta = 0.5;
OMA_PA = 0.5;

eta = 4;


pair_idx_tmp = paircombs(2*K);
pair_idx = 2*K+1-fliplr(pair_idx_tmp);
pair_idx = reshape((pair_idx)',2,K,size(pair_idx,1));
pair_idx = sort(permute(pair_idx,[2 1 3]),2);



Exhaustive_pairing = zeros(K,2,length(Pt));
RP_user_pairing = zeros(K,2,length(Pt));
User_pre_grouping = zeros(K,2,length(Pt));
User_pre_grouping_NLUPA = zeros(K,2,length(Pt));
Hungarian_pairing = zeros(K,2,length(Pt));
Simulated_Anealing_Pairing = zeros(K,2,length(Pt));
En_User_pre_grouping = zeros(K,2,length(Pt));
En_Hungarian_pairing = zeros(K,2,length(Pt));


sum_EP_opt_M_j = zeros(NNN,length(Pt));
sum_RP_opt_M_j = zeros(NNN,length(Pt));
sum_UPG_opt_M_j = zeros(NNN,length(Pt));
sum_NLUPA_opt_M_j = zeros(NNN,length(Pt));
sum_HAP_opt_M_j = zeros(NNN,length(Pt));
sum_SAP_opt_M_j = zeros(NNN,length(Pt));
sum_En_UPG_opt_M_j = zeros(NNN,length(Pt));
sum_En_HAP_opt_M_j = zeros(NNN,length(Pt));

sum_OMA_opt_M_j = zeros(NNN,length(Pt));


EP_opt_M = zeros(K,length(Pt));
RP_opt_M = zeros(K,length(Pt));
UPG_opt_M = zeros(K,length(Pt));
NULPA_opt_M = zeros(K, length(Pt));
HAP_opt_M = zeros(K, length(Pt));
SAP_opt_M = zeros(K, length(Pt));
En_UPG_opt_M = zeros(K,length(Pt));
En_HAP_opt_M = zeros(K, length(Pt));

target_BLER_j = zeros(2*K,NNN,length(Pt));

parfor u=1:length(Pt)
    for jj = 1:NNN
        h = (randn(1,N)+1i*randn(1,N));
        lamda = mean(abs(h).^2);
        % Generate user randomly
        user_distance = randi([50 500],1,2*K);
        user_distance = sort(user_distance);
        
        % Draw target BLER between 1e-8 and 1e-4 randomly
%         target_BLER = (1e-4 - 1e-8).*rand(1,2*K) + 1e-8;
        exponent = 4 + 4*rand(1, 2*K);
        target_BLER_j(:,jj,u) = 10.^-exponent;
        
%         exponent = 6 + 2*rand(1,K);
%         target_BLER_n = 10.^-exponent;
%         exponent = 4 + 2*rand(1,K);
%         target_BLER_f = 10.^-exponent;
        
%         exponent = 4 + 2*rand(1,5);
%         target_BLER_n = 10.^-exponent;
%         exponent = 6 + 2*rand(1,5);
%         target_BLER_f = 10.^-exponent;
        
%         target_BLER_j(:,jj,u) = [target_BLER_n target_BLER_f];
        target_BLER = target_BLER_j(:,jj,u);
        % Exhaustive Paring (EP)
        exhaustive_pairing = user_distance(pair_idx);
        target_BLER_EP = target_BLER(pair_idx);
        
        [sum_EP_opt_M_j(jj,u), EP_opt_M(:,u), Exhaustive_pairing(:,:,u)]=...
            EP(exhaustive_pairing, NN, K, target_BLER_EP, rho(u), eta, lamda);
        
        % Random Paring (RP)
        [sum_RP_opt_M_j(jj,u), RP_opt_M(:,u), RP_user_pairing(:,:,u)]=...
            RP(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        % User Pre-Grouping
        [sum_UPG_opt_M_j(jj,u), UPG_opt_M(:,u), User_pre_grouping(:,:,u)] =...
            UPG_opt_delta(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        % User Pre-Grouping NLUPA
        [sum_NLUPA_opt_M_j(jj,u), NULPA_opt_M(:,u), User_pre_grouping_NLUPA(:,:,u)] =...
            UPG_NLUPA(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        % Hungarian Algorithm Pairing
        [sum_HAP_opt_M_j(jj,u), HAP_opt_M(:,u), Hungarian_pairing(:,:,u)] =...
            HAP(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        % Simulated Annealing Pairing (SAP)
        [sum_SAP_opt_M_j(jj,u), SAP_opt_M(:,u), Simulated_Anealing_Pairing(:,:,u)] =...
            SAP(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        % OMA 
        [sum_OMA_opt_M_j(jj,u),~] = ...
            OMA(user_distance, NN, K, target_BLER, rho(u), beta, OMA_PA, eta, lamda);
        
        % Enhanced Hungarian Algorithm Pairing
        [sum_En_HAP_opt_M_j(jj,u), En_HAP_opt_M(:,u), En_Hungarian_pairing(:,:,u)] =...
            En_HAP(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        % Enhanced User Pre-Grouping
        [sum_En_UPG_opt_M_j(jj,u), En_UPG_opt_M(:,u), En_User_pre_grouping(:,:,u)] =...
            En_UPG_opt_delta(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
          
    end
end


sum_EP_opt_M = mean(sum_EP_opt_M_j);
sum_RP_opt_M = mean(sum_RP_opt_M_j);
sum_UPG_opt_M = mean(sum_UPG_opt_M_j);
sum_NLUPA_opt_M = mean(sum_NLUPA_opt_M_j);
sum_HAP_opt_M = mean(sum_HAP_opt_M_j);
sum_SAP_opt_M = mean(sum_SAP_opt_M_j);
sum_En_UPG_opt_M = mean(sum_En_UPG_opt_M_j);
sum_En_HAP_opt_M = mean(sum_En_HAP_opt_M_j);

sum_OMA_opt_M = mean(sum_OMA_opt_M_j);

% % Save variable
% path_str = ['C:\Users\eric7\Desktop\WeiJie\Thesis\Thesis Result\UPdata_0218'];
%  
% save(path_str, 'sum_EP_opt_M','sum_RP_opt_M', 'sum_UPG_opt_M', 'sum_NLUPA_opt_M'...
%     ,'sum_HAP_opt_M');



figure (1)

plot(Pt, sum_RP_opt_M,'--b');
hold on; grid on;
plot(Pt, sum_UPG_opt_M, 'o', 'Color',[1 0.5 0]);
plot(Pt, sum_NLUPA_opt_M, '-+m');
plot(Pt, sum_EP_opt_M, 'r');
plot(Pt, sum_HAP_opt_M, '.g');
plot(Pt, sum_SAP_opt_M, 's', 'Color', [0.3010 0.7450 0.9330]);

plot(Pt, sum_En_UPG_opt_M, '--', 'Color',[1 0.5 0]);
plot(Pt, sum_En_HAP_opt_M, '--g');

plot(Pt,sum_OMA_opt_M,'c');


xlabel('Transmitted power (dBm)');
ylabel('Blocklength (Channel use)');
legend('RP','UPG w/o Re-Grouping', 'NLUPA', ...
        'EP', 'HAP w/o Re-Grouping',...
        'SAP',...
        'UPG', 'HAP',...
        'OMA');
set(gca, 'FontName', 'Times New Roman'); 

% create a new pair of axes inside current figure
axes('position',[.65 .175 .25 .25])
box on % put box around new pair of axes
indexOfInterest = 1:3;
plot(Pt(indexOfInterest), sum_UPG_opt_M(indexOfInterest), 'o', 'Color',[1 0.5 0]);
hold on; grid on;
% plot(Pt(indexOfInterest), sum_NLUPA_opt_M(indexOfInterest), '-+m');
plot(Pt(indexOfInterest), sum_EP_opt_M(indexOfInterest), 'r');
plot(Pt(indexOfInterest), sum_HAP_opt_M(indexOfInterest), '.g');
plot(Pt(indexOfInterest), sum_SAP_opt_M(indexOfInterest), 's', 'Color', [0.3010 0.7450 0.9330]);

plot(Pt(indexOfInterest), sum_En_UPG_opt_M(indexOfInterest), '--', 'Color',[1 0.5 0]);
plot(Pt(indexOfInterest), sum_En_HAP_opt_M(indexOfInterest), '--g');

% name_str = ['UP_test_' num2str(NNN) 'times.png'];
% saveas(gcf,name_str);
% 
% saveas(gcf,'UP_test.png');