clc; clear variables; close all;

N = 1e6; % number of channel tap
NNN = 1000; % number of Monte Carlo
K = 5;  % number of cluster (number of user  = 2K)
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


EP_e_delay_j = zeros(NNN,length(Pt),2*K,2);
RP_e_delay_j = zeros(NNN,length(Pt),2*K,2);
UPG_e_delay_j = zeros(NNN,length(Pt),2*K,2);
NLUPA_e_delay_j = zeros(NNN,length(Pt),2*K,2);
HAP_e_delay_j = zeros(NNN,length(Pt),2*K,2);
SAP_e_delay_j = zeros(NNN,length(Pt),2*K,2);
En_UPG_e_delay_j = zeros(NNN,length(Pt),2*K,2);
En_HAP_e_delay_j = zeros(NNN,length(Pt),2*K,2);
OMA_e_delay_j = zeros(NNN,length(Pt),2*K,2);


EP_delay_res_1 = zeros(1,2*K,'double');
EP_delay_res_2 = zeros(1,2*K,'double');
RP_delay_res_1 = zeros(1,2*K,'double');
RP_delay_res_2 = zeros(1,2*K,'double');
UPG_delay_res_1 = zeros(1,2*K,'double');
UPG_delay_res_2 = zeros(1,2*K,'double');
NLUPA_delay_res_1 = zeros(1,2*K,'double');
NLUPA_delay_res_2 = zeros(1,2*K,'double');
HAP_delay_res_1 = zeros(1,2*K,'double');
HAP_delay_res_2 = zeros(1,2*K,'double');
SAP_delay_res_1 = zeros(1,2*K,'double');
SAP_delay_res_2 = zeros(1,2*K,'double');
En_HAP_delay_res_1 = zeros(1,2*K,'double');
En_HAP_delay_res_2 = zeros(1,2*K,'double');
En_UPG_delay_res_1 = zeros(1,2*K,'double');
En_UPG_delay_res_2 = zeros(1,2*K,'double');

sum_OMA_opt_M_j = zeros(NNN,length(Pt));


EP_opt_M_j = zeros(NNN,K,length(Pt));
RP_opt_M_j = zeros(NNN,K,length(Pt));
UPG_opt_M_j = zeros(NNN,K,length(Pt));
NULPA_opt_M_j = zeros(NNN,K, length(Pt));
HAP_opt_M_j = zeros(NNN,K, length(Pt));
SAP_opt_M_j = zeros(NNN,K, length(Pt));
En_UPG_opt_M_j = zeros(NNN,K,length(Pt));
En_HAP_opt_M_j = zeros(NNN,K, length(Pt));
OMA_opt_M_j = zeros(NNN,2*K,length(Pt));


% Near users have strict target BLER
% target_BLER = [1e-5 1e-5 1e-5 1e-5 1e-5 ...
%                1e-4 1e-4 1e-4 1e-4 1e-4];

target_BLER = [1e-7 1e-7 1e-7 1e-7 1e-7 ...
               1e-4 1e-4 1e-4 1e-4 1e-4];

% target_BLER = [1e-8 1e-8 1e-8 1e-8 1e-8 ...
%                1e-5 1e-5 1e-5 1e-5 1e-5];
           
% target_BLER = [1e-8 1e-8 1e-8 1e-8 ...
%                1e-5 1e-5 1e-5 1e-5];

% Far users have strict target BLER
% target_BLER = [1e-4 1e-4 1e-4 1e-4 1e-4 ...
%                1e-5 1e-5 1e-5 1e-5 1e-5];

% target_BLER = [1e-4 1e-4 1e-4 1e-4 1e-4 ...
%                1e-7 1e-7 1e-7 1e-7 1e-7];


for u=1:length(Pt)
    for jj = 1:NNN
        h = (randn(1,N)+1i*randn(1,N));
        lamda = mean(abs(h).^2);
        % Generate user randomly
        user_distance = randi([50 500],1,2*K);
        user_distance = sort(user_distance);
                
        % Exhaustive Paring (EP)       
        [sum_EP_opt_M_j(jj,u), EP_opt_M_j(jj,:,u), Exhaustive_pairing(:,:,u), EP_idx]=...
            EP(user_distance, pair_idx, NN, K, target_BLER, rho(u), eta, lamda);
        
        [EP_delay_tmp]=e_delay_machine(NN, 2, user_distance, EP_idx, K, rho(u), eta, target_BLER, lamda);
        EP_delay_res_1(EP_idx) = EP_delay_tmp(:,:,1);
        EP_delay_res_2(EP_idx) = EP_delay_tmp(:,:,2);
        EP_e_delay_j(jj,u,:,1) = EP_delay_res_1;
        EP_e_delay_j(jj,u,:,2) = EP_delay_res_2;
        
        % Random Paring (RP)
        [sum_RP_opt_M_j(jj,u), RP_opt_M_j(jj,:,u), RP_user_pairing(:,:,u), RP_idx]=...
            RP(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        [RP_delay_tmp]=e_delay_machine(NN, 2, user_distance, RP_idx, K, rho(u), eta, target_BLER, lamda);
        RP_delay_res_1(RP_idx) = RP_delay_tmp(:,:,1);
        RP_delay_res_2(RP_idx) = RP_delay_tmp(:,:,2);
        RP_e_delay_j(jj,u,:,1) = RP_delay_res_1;
        RP_e_delay_j(jj,u,:,2) = RP_delay_res_2;
        
%         % User Pre-Grouping
%         [sum_UPG_opt_M_j(jj,u), UPG_opt_M_j(jj,:,u), User_pre_grouping(:,:,u), UPG_idx] =...
%             UPG_opt_delta(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
%         
%         [UPG_delay_tmp]=e_delay_machine(NN, 2, user_distance, UPG_idx, K, rho(u), eta, target_BLER, lamda);
%         UPG_delay_res_1(UPG_idx) = UPG_delay_tmp(:,:,1);
%         UPG_delay_res_2(UPG_idx) = UPG_delay_tmp(:,:,2);
%         UPG_e_delay_j(jj,u,:,1) = UPG_delay_res_1;
%         UPG_e_delay_j(jj,u,:,2) = UPG_delay_res_2;      
        
        % User Pre-Grouping NLUPA
        [sum_NLUPA_opt_M_j(jj,u), NULPA_opt_M_j(jj,:,u), User_pre_grouping_NLUPA(:,:,u), NLUPA_idx] =...
            UPG_NLUPA(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        [NLUPA_delay_tmp]=e_delay_machine(NN, 2, user_distance, NLUPA_idx, K, rho(u), eta, target_BLER, lamda);
        NLUPA_delay_res_1(NLUPA_idx) = NLUPA_delay_tmp(:,:,1);
        NLUPA_delay_res_2(NLUPA_idx) = NLUPA_delay_tmp(:,:,2);
        NLUPA_e_delay_j(jj,u,:,1) = NLUPA_delay_res_1;
        NLUPA_e_delay_j(jj,u,:,2) = NLUPA_delay_res_2; 
        
%         % Hungarian Algorithm Pairing
%         [sum_HAP_opt_M_j(jj,u), HAP_opt_M_j(jj,:,u), Hungarian_pairing(:,:,u), HAP_idx] =...
%             HAP(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
%         
%         [HAP_delay_tmp]=e_delay_machine(NN, 2, user_distance, HAP_idx, K, rho(u), eta, target_BLER, lamda);
%         HAP_delay_res_1(HAP_idx) = HAP_delay_tmp(:,:,1);
%         HAP_delay_res_2(HAP_idx) = HAP_delay_tmp(:,:,2);
%         HAP_e_delay_j(jj,u,:,1) = HAP_delay_res_1;
%         HAP_e_delay_j(jj,u,:,2) = HAP_delay_res_2;
        
        % Simulated Annealing Pairing (SAP)
        [sum_SAP_opt_M_j(jj,u), SAP_opt_M_j(jj,:,u), Simulated_Anealing_Pairing(:,:,u), SAP_idx] =...
            SAP(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        [SAP_delay_tmp]=e_delay_machine(NN, 2, user_distance, SAP_idx, K, rho(u), eta, target_BLER, lamda);
        SAP_delay_res_1(SAP_idx) = SAP_delay_tmp(:,:,1);
        SAP_delay_res_2(SAP_idx) = SAP_delay_tmp(:,:,2);
        SAP_e_delay_j(jj,u,:,1) = SAP_delay_res_1;
        SAP_e_delay_j(jj,u,:,2) = SAP_delay_res_2;
        
        % OMA 
        [OMA_delay_tmp]=e_delay_machine_OMA(NN, 2, user_distance, K, rho(u), eta, beta, OMA_PA, target_BLER, lamda);
        OMA_e_delay_j(jj,u,:,1) = OMA_delay_tmp(:,1);
        OMA_e_delay_j(jj,u,:,2) = OMA_delay_tmp(:,2);
        
        % Enhanced Hungarian Algorithm Pairing
        [sum_En_HAP_opt_M_j(jj,u), En_HAP_opt_M_j(jj,:,u), En_Hungarian_pairing(:,:,u), En_HAP_idx] =...
            En_HAP(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        [En_HAP_delay_tmp]=e_delay_machine(NN, 2, user_distance, En_HAP_idx, K, rho(u), eta, target_BLER, lamda);
        En_HAP_delay_res_1(En_HAP_idx) = En_HAP_delay_tmp(:,:,1);
        En_HAP_delay_res_2(En_HAP_idx) = En_HAP_delay_tmp(:,:,2);
        En_HAP_e_delay_j(jj,u,:,1) = En_HAP_delay_res_1;
        En_HAP_e_delay_j(jj,u,:,2) = En_HAP_delay_res_2;
                
        % Enhanced User Pre-Grouping
        [sum_En_UPG_opt_M_j(jj,u), En_UPG_opt_M_j(jj,:,u), En_User_pre_grouping(:,:,u), En_UPG_idx] =...
            En_UPG_opt_delta(user_distance, NN, K, target_BLER, rho(u), eta, lamda);
        
        [En_UPG_delay_tmp]=e_delay_machine(NN, 2, user_distance, En_UPG_idx, K, rho(u), eta, target_BLER, lamda);
        En_UPG_delay_res_1(En_UPG_idx) = En_UPG_delay_tmp(:,:,1);
        En_UPG_delay_res_2(En_UPG_idx) = En_UPG_delay_tmp(:,:,2);
        En_UPG_e_delay_j(jj,u,:,1) = En_UPG_delay_res_1;
        En_UPG_e_delay_j(jj,u,:,2) = En_UPG_delay_res_2;      
    end
end

EP_e_delay = reshape(mean(EP_e_delay_j), length(Pt), 2*K, 2);
RP_e_delay = reshape(mean(RP_e_delay_j), length(Pt), 2*K, 2);
% UPG_e_delay = reshape(mean(UPG_e_delay_j), length(Pt), 2*K, 2);
% HAP_e_delay = reshape(mean(HAP_e_delay_j), length(Pt), 2*K, 2);
NLUPA_e_delay = reshape(mean(NLUPA_e_delay_j), length(Pt), 2*K, 2);
SAP_e_delay = reshape(mean(SAP_e_delay_j), length(Pt), 2*K, 2);
En_UPG_e_delay = reshape(mean(En_UPG_e_delay_j), length(Pt), 2*K, 2);
En_HAP_e_delay = reshape(mean(En_HAP_e_delay_j), length(Pt), 2*K, 2);
OMA_e_delay = reshape(mean(OMA_e_delay_j), length(Pt), 2*K, 2);

EP_e_delay_user_ave = reshape(mean(EP_e_delay, 2), length(Pt), 2);
RP_e_delay_user_ave = reshape(mean(RP_e_delay, 2), length(Pt), 2);
% UPG_e_delay_user_ave = reshape(mean(UPG_e_delay, 2), length(Pt), 2);
% HAP_e_delay_user_ave = reshape(mean(HAP_e_delay, 2), length(Pt), 2);
NLUPA_e_delay_user_ave = reshape(mean(NLUPA_e_delay, 2), length(Pt), 2);
SAP_e_delay_user_ave = reshape(mean(SAP_e_delay, 2), length(Pt),  2);
En_UPG_e_delay_user_ave = reshape(mean(En_UPG_e_delay, 2), length(Pt), 2);
En_HAP_e_delay_user_ave = reshape(mean(En_HAP_e_delay, 2), length(Pt), 2);
OMA_e_delay_user_ave = reshape(mean(OMA_e_delay, 2), length(Pt), 2);

ret1 = [EP_e_delay_user_ave(:,1) EP_e_delay_user_ave(:,2) ...
       RP_e_delay_user_ave(:,1) RP_e_delay_user_ave(:,2) ...
       NLUPA_e_delay_user_ave(:,1) NLUPA_e_delay_user_ave(:,2) ...
       SAP_e_delay_user_ave(:,1) SAP_e_delay_user_ave(:,2) ...
       En_UPG_e_delay_user_ave(:,1) En_UPG_e_delay_user_ave(:,2) ...
       En_HAP_e_delay_user_ave(:,1) En_HAP_e_delay_user_ave(:,2) ...
       OMA_e_delay_user_ave(:,1) OMA_e_delay_user_ave(:,2)];

 figure (1)
 b = bar(ret1); 
 b(1).FaceColor = 'r';
 b(2).FaceColor = 'k';
 b(3).FaceColor = 'b';
 b(4).FaceColor = [0.9290 0.6940 0.1250];
 b(5).FaceColor = 'm';
 b(6).FaceColor = [0.4660 0.6740 0.1880];
 b(7).FaceColor = [0.3010 0.7450 0.9330];
 b(8).FaceColor = [0.6350 0.0780 0.1840];
 b(9).FaceColor = [1 0.5 0];
 b(10).FaceColor = [0 0.4470 0.7410];
 b(11).FaceColor = 'g';
 b(12).FaceColor = [0.8500 0.3250 0.0980];
 b(13).FaceColor = 'c';
 b(14).FaceColor = [0.4940 0.1840 0.5560];
 legend('EP 1 shot', 'EP HARQ(2)',...
        'RP 1 shot', 'RP HARQ(2)',...
        'NLUPA 1 shot', 'NLUPA HARQ(2)',...
        'SAP 1 shot', 'SAP HARQ(2)',...
        'UPG 1 shot', 'UPG HARQ(2)',...
        'HAP 1 shot', 'HAP HARQ(2)',...
        'OMA 1 shot', 'OMA HARQ(2)');
    
xlabel('Transmitted power level');
ylabel('Expected delay (Channel use)');
set(gca, 'FontName', 'Times New Roman'); 
    
figure (2) 
ret2 = [EP_e_delay(6,:,1)' RP_e_delay(6,:,1)'...
     NLUPA_e_delay(6,:,1)' SAP_e_delay(6,:,1)'...
     En_UPG_e_delay(6,:,1)' En_HAP_e_delay(6,:,1)'...
     OMA_e_delay(6,:,1)'];
 bb = bar(ret2);
 bb(1).FaceColor = 'r';
 bb(2).FaceColor = 'b';
 bb(3).FaceColor = 'm';
 bb(4).FaceColor = [0.3010 0.7450 0.9330];
 bb(5).FaceColor = [1 0.5 0];
 bb(6).FaceColor = 'g';
 bb(7).FaceColor = 'c';
 
legend('EP 1 shot',...
    'RP 1 shot',...
    'NLUPA 1 shot',...
    'SAP 1 shot',...
    'UPG 1 shot',...
    'HAP 1 shot',...
    'OMA 1 shot');

xlabel('User');
ylabel('Expected delay (Channel use)');
set(gca, 'FontName', 'Times New Roman'); 
 
figure (3) 
ret3 = [EP_e_delay(6,:,2)' RP_e_delay(6,:,2)'...
     NLUPA_e_delay(6,:,2)' SAP_e_delay(6,:,2)'...
     En_UPG_e_delay(6,:,2)' En_HAP_e_delay(6,:,2)'...
     OMA_e_delay(6,:,2)'];
 bbb = bar(ret3);
 bbb(1).FaceColor = 'r';
 bbb(2).FaceColor = 'b';
 bbb(3).FaceColor = 'm';
 bbb(4).FaceColor = [0.3010 0.7450 0.9330];
 bbb(5).FaceColor = [1 0.5 0];
 bbb(6).FaceColor = 'g';
 bbb(7).FaceColor = 'c';
 
 legend('EP HARQ(2)',...
        'RP HARQ(2)',...
        'NLUPA HARQ(2)',...
        'SAP HARQ(2)',...
        'UPG HARQ(2)',...
        'HAP HARQ(2)',...
        'OMA HARQ(2)');
    
xlabel('User');
ylabel('Expected delay (Channel use)');
set(gca, 'FontName', 'Times New Roman'); 
 
 
sum_EP_opt_M = mean(sum_EP_opt_M_j);
sum_RP_opt_M = mean(sum_RP_opt_M_j);
% sum_UPG_opt_M = mean(sum_UPG_opt_M_j);
sum_NLUPA_opt_M = mean(sum_NLUPA_opt_M_j);
% sum_HAP_opt_M = mean(sum_HAP_opt_M_j);
sum_SAP_opt_M = mean(sum_SAP_opt_M_j);
sum_En_UPG_opt_M = mean(sum_En_UPG_opt_M_j);
sum_En_HAP_opt_M = mean(sum_En_HAP_opt_M_j);

