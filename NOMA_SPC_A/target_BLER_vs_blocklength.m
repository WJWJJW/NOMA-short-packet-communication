clc; clear variables; close all;

N = 1e6; % number of channel tap
NNN = 1000; % number of Monte Carlo
K = 5  % number of cluster (number of user  = 2K)
NN = 256; % number of information bit
N1 = NN;
N2 = NN;


Pt = 30;               %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

% AWGN
% BW = 10^7;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt./ no;
RHO = pow2db(rho);

eta = 4;

target_BLER_i = [1e-8 1e-7 1e-6 1e-5 1e-4 1e-3];
target_BLER_j = [1e-8 1e-7 1e-6 1e-5 1e-4 1e-3];

user_distance = [100 200];

sum_opt_M1 = zeros(1,length(target_BLER_i));
sum_opt_M2 = zeros(1,length(target_BLER_i));
cur_opt_M1 = zeros(2,length(target_BLER_i));
cur_opt_M2 = zeros(2,length(target_BLER_i));

for u=1:length(target_BLER_i)
    h = (randn(1,N)+1i*randn(1,N));
    lamda = mean(abs(h).^2);
    [sum_opt_M1(u), cur_opt_M1(:,u)] = ...
        M_cal_Mod(N1,user_distance,1,[target_BLER_i(1),target_BLER_j(u)],rho,eta,lamda);
    [sum_opt_M2(u), cur_opt_M2(:,u)] = ...
        M_cal_Mod(N1,user_distance,1,[target_BLER_i(u),target_BLER_j(1)],rho,eta,lamda);
end



figure(1)
plot(1:length(target_BLER_i), sum_opt_M1,'r');
hold on ; grid on;
plot(1:length(target_BLER_i), sum_opt_M2,'b');

