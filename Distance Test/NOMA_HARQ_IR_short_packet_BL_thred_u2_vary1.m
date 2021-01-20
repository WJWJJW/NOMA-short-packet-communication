% This script analysis the relation between blocklength , power allocation
% coefficient and user 2 distance
% calculate blocklength and power allocation coefficient according
% different conditon

clc; clear variables; close all;
N = 1e6;
N1 = 80;
N2 = 80;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = -4;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt/ no;

eta = 4;

d1 = 120;
pair_d2_UPG = 188;
pair_d2_EP = 167;
% diff = 10:2:300;
% d2 = d1+diff;
d2 = 150:2:330;
check1 = zeros(length(d2),1);
check2 = zeros(length(d2),1);

delta = 1/1.2;

h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
lamda1 = mean(abs(h1).^2);

for dd=1:length(d2)
    h2 = sqrt(1/2*d2(dd)^-eta)*(randn(1,N)+1i*randn(1,N));
    lamda2 = mean(abs(h2).^2);
    
    contraint(dd) = (lamda1*eplsion1R)/(lamda2*eplsion2R);
    k = contraint(dd)*delta;
   
    % optimal alpha1 and blocklength calculation
    if d2(dd)/d1 > (eplsion2R/eplsion1R)^0.25 && (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - (lamda2*eplsion2R.*rho + 4)*lamda2*eplsion2R >0
        opt_a1(dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
                  4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R))) ...
                  / (2*lamda2*eplsion2R*rho*(lamda1*eplsion1R-lamda2*eplsion2R));
        opt_M(dd) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R)))...
                             / (2*lamda2*eplsion2R)));
    
    elseif d2(dd)/d1 > (eplsion2R/eplsion1R)^0.25 && (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - (lamda2*eplsion2R.*rho + 4)*lamda2*eplsion2R <0
        opt_a1(dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
                  4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R))) ...
                  / (2*k*lamda2*eplsion2R*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R));
        opt_M(dd) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R)))...
                             / (2*k*lamda2*eplsion2R)));
        check1(dd) = 1;
    else
        opt_a1(dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
                  4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R))) ...
                  / (2*k*lamda2*eplsion2R*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R));
        opt_M(dd) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R)))...
                             / (2*k*lamda2*eplsion2R)));
        check2(dd) = 1;
    end

end

figure (1)
syms D2
eqn_constraint1 = d1^-eta*eplsion1R - D2^-eta*eplsion2R == 0;
D2_sol1 = vpasolve(eqn_constraint1, D2);
eqn_constraint2 = (D2^-eta*eplsion2R*rho + 2)*d1^-eta*eplsion1R - (D2^-eta.*eplsion2R.*rho + 4)*D2^-eta*eplsion2R == 0;
D2_sol2 = vpasolve(eqn_constraint2, D2);

yyaxis left
plot(d2, real(opt_a1),'b');
hold on; grid on;

yline(0.5,'-');

ylabel('\alpha_1');
axis([-inf inf 0 1]);

yyaxis right
plot(d2, real(opt_M), 'o', 'Color',[1 0.5 0]);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('blocklength');


legend('opt \alpha_1','blocklength');
    
xline(double(D2_sol1(2)),'-','Threshold1');
xline(double(D2_sol2(2)),'-.','Threshold2');
xline(pair_d2_UPG,'--','PairUPG d_2');
xline(pair_d2_EP,'--','PairEP d_2');
% axis([-inf inf 0 2000]);


xlabel('user 2 distance');



