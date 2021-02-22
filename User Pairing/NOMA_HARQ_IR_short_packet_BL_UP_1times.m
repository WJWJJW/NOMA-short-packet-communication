%  This script analysis performance of UPG EP RP and HAP 1 time
clc; clear variables; close all;


N = 1e6; % number of Monte Carlo
K = 5;  % number of cluster (number of user  = 2K)
NN = 256; % number of information bit
N1 = NN;
N2 = NN;

eplsion1R = 10^-5;
eplsion2R = 10^-4;


Pt = 20:2:30;                    %Transmit Power in dBm
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

Exhaustive_pairing = zeros(K,2,length(Pt));
RP_user_pairing = zeros(K,2,length(Pt));
User_pre_grouping = zeros(K,2,length(Pt));
User_pre_grouping_NLUPA = zeros(K,2,length(Pt));
Hungarian_pairing = zeros(K,2,length(Pt));

sum_EP_opt_M = zeros(1,length(Pt));
sum_RP_opt_M = zeros(1,length(Pt));
sum_UPG_opt_M = zeros(1,length(Pt));
sum_NLUPA_opt_M = zeros(1,length(Pt));
sum_HAP_opt_M = zeros(1,length(Pt));

sum_OMA_opt_M_a = zeros(1,length(Pt));
sum_OMA_opt_M_b = zeros(1,length(Pt));
sum_OMA_opt_M_c = zeros(1,length(Pt));
sum_OMA_opt_M_d = zeros(1,length(Pt));

EP_opt_M = zeros(K,length(Pt));
RP_opt_M = zeros(K,length(Pt));
UPG_opt_M = zeros(K,length(Pt));
NULPA_opt_M = zeros(K, length(Pt));
HAP_opt_M = zeros(K, length(Pt));

OMA_opt_M_a = zeros(2*K, length(Pt));
OMA_opt_M_b = zeros(2*K, length(Pt));
OMA_opt_M_c = zeros(2*K, length(Pt));
OMA_opt_M_d = zeros(2*K, length(Pt));

dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);

user_distance = randi([50 300],1,2*K);
user_distance = sort(user_distance);


pair_idx_tmp = paircombs(2*K);
pair_idx = 2*K+1-fliplr(pair_idx_tmp);
clear pair_idx_tmp;
exhaustive_pairing = reshape(user_distance(pair_idx)',K,2,length(pair_idx));
clear pair_idx;



for u=1:length(Pt)
    h = (randn(1,N)+1i*randn(1,N));
    lamda = mean(abs(h).^2);
    % Exhaustive Paring (EP)
    [sum_EP_opt_M(u), EP_opt_M(:,u), Exhaustive_pairing(:,:,u)]=...
        EP(exhaustive_pairing, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda);
    % Random Paring (RP)
    [sum_RP_opt_M(u), RP_opt_M(:,u), RP_user_pairing(:,:,u)]=...
        RP(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda);
    % User Pre-Grouping
    [sum_UPG_opt_M(u), UPG_opt_M(:,u), User_pre_grouping(:,:,u)] =...
        UPG_opt_delta(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda);
%     % Simulated Annealing Pairing
%     [sum_SAP_opt_M(u), SAP_opt_M(:,u), Simulated_Anealing_Pairing(:,:,u)] =...
%         SAP(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda, delta);
    
    % Hungarian Algorithm Pairing
    [sum_HAP_opt_M(u), HAP_opt_M(:,u), Hungarian_pairing(:,:,u)] =...
        HAP(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda);
    
    % User Pre-Grouping NLUPA
    [sum_NLUPA_opt_M(u), NULPA_opt_M(:,u), User_pre_grouping_NLUPA(:,:,u)] =...
        UPG_NLUPA(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda);
    
    % OMA (equal reliability constraint)
    [sum_OMA_opt_M_a(u), OMA_opt_M_a(:,u)] = ...
        OMA_a(user_distance, NN, K, eplsion1R, eplsion1R, rho(u), beta, OMA_PA, eta, lamda);
    [sum_OMA_opt_M_b(u), OMA_opt_M_b(:,u)] = ...
        OMA_a(user_distance, NN, K, eplsion2R, eplsion2R, rho(u), beta, OMA_PA, eta, lamda);
    % OMA (non-equal reliability constraint)
    [sum_OMA_opt_M_c(u), OMA_opt_M_c(:,u)] = ...
        OMA_a(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), beta, OMA_PA, eta, lamda);
    % OMA (reliability constraint according to pair rule)
    [sum_OMA_opt_M_d(u), OMA_opt_M_d(:,u)] = ...
        OMA_b(User_pre_grouping(:,:,u), NN, K, eplsion1R, eplsion2R, rho(u), beta, OMA_PA, eta, lamda);
end


% % Save variable
% test_idx = 1;
% path_str1 = ['E:\WeiJie\NOMA\Matlab\Thesis_data\0105_UP_test\UP_1time_' num2str(test_idx) '_u_distribution'];
% path_str2 = ['E:\WeiJie\NOMA\Matlab\Thesis_data\0105_UP_test\UP_1time_' num2str(test_idx) '_check'];
% path_str3 = ['E:\WeiJie\NOMA\Matlab\Thesis_data\0105_UP_test\UP_1time_' num2str(test_idx) '_u_combination'];
% 
% save(path_str1,'user_distance');
% save(path_str2, 'SAP_thred1_check','SAP_thred2_check', 'UPG_thred1_check', 'UPG_thred2_check');
% save(path_str3,'Simulated_Anealing_Pairing','User_pre_grouping');
% name_str = ['E:\WeiJie\NOMA\Matlab\Thesis_data\0105_UP_test\UP_test_' num2str(test_idx) '.png'];


figure (1)

plot(Pt, sum_RP_opt_M,'b');
hold on; grid on;
plot(Pt, sum_UPG_opt_M, 'Color',[1 0.5 0]);
plot(Pt, sum_NLUPA_opt_M, 'm');
plot(Pt, sum_EP_opt_M, 'ro');
plot(Pt, sum_HAP_opt_M, 'g*');

plot(Pt,sum_OMA_opt_M_a,'c');
plot(Pt,sum_OMA_opt_M_b,'*c');
plot(Pt,sum_OMA_opt_M_c,'oc');


ylabel('blocklength');
legend('Random Pairing','User Pre-Grouping', 'User Pre-Grouping NLUPA', ...
        'Exhaustive Paring', 'Hungarian Pairing',...
        'OMA equal relibility constraint (high)', ...
        'OMA equal relibility constraint (low)',...
        'OMA non-equal relibility constraint');
set(gca, 'FontName', 'Times New Roman'); 

figure (2)

plot(Pt, sum_RP_opt_M,'b');
hold on; grid on;
plot(Pt, sum_UPG_opt_M, 'Color',[1 0.5 0]);
plot(Pt, sum_NLUPA_opt_M, 'm');
plot(Pt, sum_EP_opt_M, 'ro');
plot(Pt, sum_HAP_opt_M, 'g*');

plot(Pt,sum_OMA_opt_M_c,'c');
plot(Pt,sum_OMA_opt_M_d,'*c');


ylabel('blocklength');
legend('Random Pairing','User Pre-Grouping', 'User Pre-Grouping NLUPA', ...
        'Exhaustive Paring', 'Hungarian Pairing',...
        'OMA non-equal relibility constraint',...
        'OMA Pair');

set(gca, 'FontName', 'Times New Roman');

% saveas(gcf,name_str);