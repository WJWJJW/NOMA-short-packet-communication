clc; clear variables; close all;
N = 1e6;
eta = 4;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = [20 30 40];                    %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rrho = pt./ no;

d2 = 10:100:25000;
syms D1
for rr = 1:length(rrho)
    rho = rrho(rr);
    for dd = 1:length(d2)
        eqn_constraint1 = D1^-eta*eplsion1R - d2(dd)^-eta*eplsion2R == 0;
        D1_sol1(rr,dd,:) = vpasolve(eqn_constraint1, D1);
        eqn_constraint2 = (d2(dd)^-eta*eplsion2R*rho + 2)*D1^-eta*eplsion1R - (d2(dd)^-eta.*eplsion2R.*rho + 4)*d2(dd)^-eta*eplsion2R == 0;
        D1_sol2(rr,dd,:) = vpasolve(eqn_constraint2, D1);
    end

end






figure (1)
yyaxis left
plot(d2, D1_sol2(1,:,2),'b', 'linewidth', 1.5);
hold on; grid on;

plot(d2, D1_sol2(2,:,2),'ob');

plot(d2, D1_sol2(3,:,2),'*b');
ylabel('d1');


yyaxis right
plot(d2, d2./D1_sol2(1,:,2),'Color',[1 0.5 0], 'linewidth', 1.5);
plot(d2, d2./D1_sol2(2,:,2),'c', 'linewidth', 1.5);
plot(d2, d2./D1_sol2(3,:,2),'k', 'linewidth', 1.5);
ylabel('ratio of distance');

legend('20 dBm','30 dBm','40 dBm',...
       'ratio for 20 dBm','ratio for 30 dBm','ratio for 40 dBm');
xlabel('d2');


    