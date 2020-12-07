clc; clear variables; close all;
N = 1e6;
eta = 4;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 40;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt/ no;

d1 = 1:1:500;
syms D2
for dd = 1:length(d1)
    eqn_constraint1 = d1(dd)^-eta*eplsion1R - D2^-eta*eplsion2R == 0;
    D2_sol1(dd,:) = vpasolve(eqn_constraint1, D2);
    eqn_constraint2 = (D2^-eta*eplsion2R*rho + 2)*d1(dd)^-eta*eplsion1R - (D2^-eta.*eplsion2R.*rho + 4)*D2^-eta*eplsion2R == 0;
    D2_sol2(dd,:) = vpasolve(eqn_constraint2, D2);
end

figure (1)
yyaxis left
plot(d1, D2_sol1(:,2),'r', 'linewidth', 1.5);
hold on; grid on;
plot(d1, D2_sol2(:,2),'b', 'linewidth', 1.5);
ylabel('d2');


yyaxis right

plot(d1, D2_sol2(:,2)./d1.');
ylabel('ratio of distance');

% legend('constraint 1','constraint 2','difference of distance');
xlabel('d1');
