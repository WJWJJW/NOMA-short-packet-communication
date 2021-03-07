% This script analysis the relation between blocklength , power allocation
% coefficient and user 2 distance
% calculate blocklength and power allocation coefficient according
% different conditon

clc; clear variables; close all;
N = 1e7;
N1 = 256;
N2 = 256;

eplsion1R = 10^-4;
eplsion2R = 10^-5;

% Tx power
Pt = 30;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

% AWGN
% BW = 10^6;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt/ no;

eta = 4;

d1 = 100;
% pair_d2_UPG = 188;
% pair_d2_EP = 167;
% diff = 10:2:300;
% d2 = d1+diff;
d2 = 120:2:300;


h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
lamda1 = mean(abs(h1).^2);

% delta = 1/2;
[delta] = delta_finder (lamda1*rho*eplsion1R);

for dd=1:length(d2)
    h2 = sqrt(1/2*d2(dd)^-eta)*(randn(1,N)+1i*randn(1,N));
    lamda2 = mean(abs(h2).^2);
    
    contraint = (lamda1*eplsion1R)/(lamda2*eplsion2R);
    theta = contraint * delta;
    
    % Adjust target BLER according to condition 2
    check = (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - (lamda2*eplsion2R.*rho + 4)*lamda2*eplsion2R > 0;
    
    switch (check)
        case 0
            eplsion2R_m = theta*eplsion2R;
        case 1
            eplsion2R_m = eplsion2R;
    end
    
    % optimal alpha1 and blocklength calculation
    opt_a1(dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
              4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m))) ...
              / (2*lamda2*eplsion2R_m*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m));
    opt_M(dd) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)))...
                         / (2*lamda2*eplsion2R_m)));
end

figure (1)
syms D2
% eqn_constraint1 = d1^-eta*eplsion1R - D2^-eta*eplsion2R == 0;
% D2_sol1 = vpasolve(eqn_constraint1, D2);
eqn_constraint2 = (D2^-eta*eplsion2R*rho + 2)*d1^-eta*eplsion1R - (D2^-eta.*eplsion2R.*rho + 4)*D2^-eta*eplsion2R == 0;
D2_sol2 = vpasolve(eqn_constraint2, D2);

yyaxis left
plot(d2, real(opt_a1),'b');
hold on; grid on;

yline(0.5,'-');

ylabel('\alpha_1');
axis([-inf inf 0 0.1]);

yyaxis right
plot(d2, real(opt_M), 'Color',[1 0.5 0]);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('Blocklength (Channel uses)');


legend('\alpha_1^{opt}','M_{opt}^{1 pair}');
    
% xline(double(D2_sol1(2)),'-','Threshold1');
% % xx = xline(double(D2_sol2(2)),'-.','Bound (38)','HandleVisibility','off');
% xx.FontName = 'Times New Roman';
% xline(pair_d2_UPG,'--','PairUPG d_2');
% xline(pair_d2_EP,'--','PairEP d_2');
% axis([-inf inf 0 2000]);


xlabel('Distance of far user (meter)');

set(gca, 'FontName', 'Times New Roman');

