clc; clear variables; close all;
% Transmitted Power
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

beta = 0.5;
OMA_PA = 0.5;

eta = 4;
target_BLER_n = [1e-8 1e-6];
target_BLER_f = 1e-8:1e-6:1e-4;

user_distance = [100 200];

d_min = zeros(length(target_BLER_n),length(target_BLER_f),length(user_distance));

for n=1:length(target_BLER_n)
    for f=1:length(target_BLER_f)
        for d=1:length(user_distance)
        d_min(n,f,d) = min_d_finder (user_distance(d),target_BLER_n(n),target_BLER_f(f),rho);
        end
    end
end

figure(1)
semilogx(target_BLER_f,d_min(1,:,1),'b', 'linewidth', 1.5);
hold on;grid on;
semilogx(target_BLER_f,d_min(1,:,2),'--b', 'linewidth', 1.5);
semilogx(target_BLER_f,d_min(2,:,1),'r', 'linewidth', 1.5);
semilogx(target_BLER_f,d_min(2,:,2),'--r', 'linewidth', 1.5);

x1 = xlabel('Target BLER for FU, $\tilde{\epsilon_{j}}$');
set(x1,'interpreter','Latex','FontSize',14);

y1 = ylabel('Min. distance for FU');
set(y1,'FontSize',14);

h = legend('$\tilde{\epsilon_{i}} = 10^{-8}, d_i = 100$',...
       '$\tilde{\epsilon_{i}} = 10^{-8}, d_i = 200$',...
       '$\tilde{\epsilon_{i}} = 10^{-6}, d_i = 100$',...
       '$\tilde{\epsilon_{i}} = 10^{-6}, d_i = 200$');
set(h,'interpreter','Latex','FontSize',14);

set(gca, 'FontName', 'Times New Roman'); 

figure(2)
semilogx(target_BLER_f,d_min(1,:,1)-user_distance(1),'b', 'linewidth', 1.5);
hold on;grid on;
semilogx(target_BLER_f,d_min(1,:,2)-user_distance(2),'--b', 'linewidth', 1.5);
semilogx(target_BLER_f,d_min(2,:,1)-user_distance(1),'r', 'linewidth', 1.5);
semilogx(target_BLER_f,d_min(2,:,2)-user_distance(2),'--r', 'linewidth', 1.5);

xlabel('Target BLER for FU');
ylabel('Difference between Min. distance for FU and distance of NU');
h2 = legend('$\tilde{\epsilon_{i}} = 10^{-8}, d_i = 100$',...
       '$\tilde{\epsilon_{i}} = 10^{-8}, d_i = 200$',...
       '$\tilde{\epsilon_{i}} = 10^{-6}, d_i = 100$',...
       '$\tilde{\epsilon_{i}} = 10^{-6}, d_i = 200$');
h2.Interpreter = 'latex';
set(gca, 'FontName', 'Times New Roman'); 


