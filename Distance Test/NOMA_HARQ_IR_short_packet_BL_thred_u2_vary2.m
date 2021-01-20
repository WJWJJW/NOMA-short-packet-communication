% This script analysis the relation between blocklength , power allocation
% coefficient and user 2 distance
% calculate blocklength and power allocation coefficient in same condition
% (same delta value)

clc; clear variables; close all;
N = 1e6;
N1 = 256;
N2 = 256;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 20;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt/ no;

eta = 4;

d1 = 50;
diff = 10:2:200;
d2 = d1+diff;


delta = ((d2(end)^-eta)*eplsion2R)/((d1^-eta)*eplsion1R);

h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
lamda1 = mean(abs(h1).^2);

for dd=1:length(diff)
    h2 = sqrt(1/2*d2(dd)^-eta)*(randn(1,N)+1i*randn(1,N));
    lamda2 = mean(abs(h2).^2);
    
    contraint(dd) = (lamda1*eplsion1R)/(lamda2*eplsion2R);
    theta = contraint(dd)*delta;
   
    opt_a1(dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
              4*(theta*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-theta*lamda2*eplsion2R))) ...
              / (2*theta*lamda2*eplsion2R*rho*(lamda1*eplsion1R-theta*lamda2*eplsion2R));
    opt_M(dd) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(theta*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-theta*lamda2*eplsion2R)))...
                         / (2*theta*lamda2*eplsion2R)));

%     opt_a1(dd) = (-1+sqrt(1+4*delta^2*rho*lamda1*eplsion1R*(1-delta)))/(2*delta*rho*lamda1*eplsion1R*(1-delta));
%     opt_M(dd) = N1/(log2(1+(-1+sqrt(1+4*delta^2*rho*lamda1*eplsion1R*(1-delta)))/(2*delta)));
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
% axis([-inf inf 0 1]);

yyaxis right
plot(d2, real(opt_M), 'o', 'Color',[1 0.5 0]);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('blocklength');

legend('opt \alpha_1','blocklength');
    
xline(double(D2_sol1(2)),'-','Threshold1');
xline(double(D2_sol2(2)),'-.','Threshold2');
% axis([-inf inf 0 2000]);


xlabel('user 2 distance');



