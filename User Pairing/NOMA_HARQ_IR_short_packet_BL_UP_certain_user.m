% This script analysis performance of EP, RP, UPG, SAP
% Certain user distance

clc; clear variables; close all;

N = 1e6; % number of Monte Carlo
K = 5;  % number of cluster (number of user  = 2K)
NN = 80; % number of information bit
N1 = NN;
N2 = NN;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 0:2:20;                    %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt./ no;
RHO = pow2db(rho);

eta = 4;

Exhaustive_pairing = zeros(K,2,length(Pt));
RP_user_pairing = zeros(K,2,length(Pt));
User_pre_grouping = zeros(K,2,length(Pt));
Simulated_Anealing_Pairing = zeros(K,2,length(Pt));

sum_EP_opt_M = zeros(1,length(Pt));
sum_RP_opt_M = zeros(1,length(Pt));
sum_UPG_opt_M = zeros(1,length(Pt));
sum_SAP_opt_M = zeros(1,length(Pt));

EP_opt_M = zeros(K,length(Pt));
RP_opt_M = zeros(K,length(Pt));
UPG_opt_M = zeros(K,length(Pt));
SAP_opt_M = zeros(K, length(Pt));

dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);

% load_idx = 2;
% load_str1 = ['E:\WeiJie\NOMA\Matlab\Thesis_data\1217_UP_test\UP_1time_' num2str(load_idx) '_u_distribution'];
% load(load_str1, 'user_distance');

% user_distance = [10 20 30 40 50 250 260 270 280 290];
% user_distance = [10 20 130 140 150 160 170 180 190 290];
% user_distance = [64 87 90 108 119 154 167 188 266 328];
% user_distance = [10 11 12 13 14 15 16 17 18 19];
user_distance = [50 51 52 53 54 55 56 57 58 59];

pair_idx_tmp = paircombs(2*K);
pair_idx = 2*K+1-fliplr(pair_idx_tmp);
clear pair_idx_tmp;
exhaustive_pairing = reshape(user_distance(pair_idx)',K,2,length(pair_idx));
clear pair_idx;

delta = 1/2;

for u=1:length(Pt)
    h = (randn(1,N)+1i*randn(1,N));
    lamda = mean(abs(h).^2);
    % Exhaustive Paring (EP)
    [sum_EP_opt_M(u), EP_opt_M(:,u), Exhaustive_pairing(:,:,u)]=...
        EP(exhaustive_pairing, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda, delta);
    % Random Paring (RP)
    [sum_RP_opt_M(u), RP_opt_M(:,u), RP_user_pairing(:,:,u)]=...
        RP(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda, delta);
    % User Pre-Grouping
    [sum_UPG_opt_M(u), UPG_opt_M(:,u), User_pre_grouping(:,:,u)] =...
        UPG(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda, delta);
    % Simulated Annealing Pairing
    [sum_SAP_opt_M(u), SAP_opt_M(:,u), Simulated_Anealing_Pairing(:,:,u)] =...
        SAP(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda, delta);
end




% % Save variable
% test_idx = 5;
% path_str1 = ['E:\WeiJie\NOMA\Matlab\Thesis_data\1217_UP_test\UP_1time_' num2str(test_idx) '_u_distribution'];
% path_str2 = ['E:\WeiJie\NOMA\Matlab\Thesis_data\1217_UP_test\UP_1time_' num2str(test_idx) '_check'];
% path_str3 = ['E:\WeiJie\NOMA\Matlab\Thesis_data\1217_UP_test\UP_1time_' num2str(test_idx) '_u_combination'];
% 
% save(path_str1,'user_distance');
% save(path_str2, 'SAP_thred1_check','SAP_thred2_check', 'UPG_thred1_check', 'UPG_thred2_check');
% save(path_str3,'Simulated_Anealing_Pairing','User_pre_grouping');
% name_str = ['UP_test_' num2str(test_idx) '.png'];


figure (1)

plot(Pt, sum_RP_opt_M,'b');
hold on; grid on;
plot(Pt, sum_UPG_opt_M,'-.cs');
plot(Pt, sum_SAP_opt_M,'Color',[1 0.5 0]);
plot(Pt, sum_EP_opt_M, 'm');

ylabel('blocklength');
legend('Random Pairing','User Pre-Grouping', 'Simulated Anealing Based Pairing', 'Exhaustive Paring');

figure (2)

plot(Pt, sum_UPG_opt_M,'-.cs');
hold on; grid on;
plot(Pt, sum_EP_opt_M, 'm');

ylabel('blocklength');
legend('User Pre-Grouping','Exhaustive Paring');
   
% saveas(gcf,name_str);