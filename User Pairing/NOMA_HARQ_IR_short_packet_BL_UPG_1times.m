% This script analysis performance of UPG EP and NLUPA

clc; clear variables; close all;
N = 1e6; % number of Monte Carlo
K = 5;  % number of cluster (number of user  = 2K)
NN = 80; % number of information bit
N1 = NN;
N2 = NN;

eplsion1R = 10^-5;
eplsion2R = 10^-4;


Pt = 10:2:20;                    %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt./ no;
RHO = pow2db(rho);

eta = 4;
dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);

user_distance = randi([10 330],1,2*K);
user_distance = sort(user_distance);


pair_idx_tmp = paircombs(2*K);
pair_idx = 2*K+1-fliplr(pair_idx_tmp);
clear pair_idx_tmp;
exhaustive_pairing = reshape(user_distance(pair_idx)',K,2,length(pair_idx));
clear pair_idx;


Exhaustive_pairing = zeros(K,2,length(Pt));
User_pre_grouping = zeros(K,2,length(Pt));
User_pre_grouping_NLUPA = zeros(K,2,length(Pt));



sum_EP_opt_M = zeros(1,length(Pt));
sum_UPG_opt_M = zeros(1,length(Pt));
sum_NLUPA_opt_M = zeros(1,length(Pt));


EP_opt_M = zeros(K,length(Pt));
UPG_opt_M = zeros(K,length(Pt));
NULPA_opt_M = zeros(K,length(Pt));


delta = 1/2;
for u=1:length(Pt)
    h = (randn(1,N)+1i*randn(1,N));
    lamda = mean(abs(h).^2);
    % Exhaustive Paring (EP)
    [sum_EP_opt_M(u), EP_opt_M(:,u), Exhaustive_pairing(:,:,u)]=...
        EP(exhaustive_pairing, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda, delta);
    % User Pre-Grouping
    [sum_UPG_opt_M(u), UPG_opt_M(:,u), User_pre_grouping(:,:,u)] =...
        UPG(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda, delta);
    % User Pre-Grouping NLUPA
    [sum_NLUPA_opt_M(u), NULPA_opt_M(:,u), User_pre_grouping_NLUPA(:,:,u)] =...
        UPG_NLUPA(user_distance, NN, K, eplsion1R, eplsion2R, rho(u), eta, lamda, delta);
end



figure (1)

plot(Pt, sum_UPG_opt_M,'b');
hold on; grid on;
plot(Pt, sum_NLUPA_opt_M,'r');
plot(Pt, sum_EP_opt_M, 'm');

ylabel('blocklength');
legend('User Pre-Grouping (Proposed)','User Pre-Grouping (NLUPA)','Exhaustive Paring');





   