clc; clear variables; close all;
N = 1e6;
N1 = 80;
N2 = 80;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 20;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt/ no;

eta = 4;
dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);
dis_thred2 = (1.64*eplsion2R/eplsion1R)^(1/eta);
dis_thred3 = (2*eplsion2R/eplsion1R)^(1/eta);


d1 = 119;
diff = 10:2:150;
d2 = d1+diff;

syms otp_a1
check1 = zeros(length(diff),1);
check2 = zeros(length(diff),1);

for dd=1:length(diff)
    
    h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
    h2 = sqrt(1/2*d2(dd)^-eta)*(randn(1,N)+1i*randn(1,N));

    lamda1 = mean(abs(h1).^2);
    lamda2 = mean(abs(h2).^2);
    
    contraint(dd) = (lamda1*eplsion1R)/(lamda2*eplsion2R);
    k = contraint(dd)/6;
    
    % solve by matlab
%     xo2 = (1+eplsion2R*lamda2*rho)/(1+eplsion2R*lamda2*rho*otp_a1);
%     xo1 = 1+otp_a1*rho*(lamda1*eplsion1R-lamda2*eplsion2R);
%     eqn_opt = (xo2) == (xo1);
%     opt_a1_mat(dd,:) = vpasolve(eqn_opt,otp_a1);

    
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

% figure (1)
% ratio = d2./d1;
% yyaxis left
% plot(ratio, real(opt_a1),'ob');
% hold on; grid on;
% % plot(ratio, real(opt_a1_mat(:,2)),'b','linewidth', 1.5);
% % plot(ratio, imag(opt_a1_mat(:,2)),'r','linewidth', 1.5);
% yline(0.5,'-');
% 
% ylabel('\alpha_1');
% axis([-inf inf 0 1]);
% 
% yyaxis right
% plot(ratio, real(opt_M),'Color',[1 0.5 0]);
% hold on; grid on;
% % plot(ratio, real(opt_M_u1),'o','Color',[1 0.5 0]);
% % plot(ratio, real(opt_M_u2),'+','Color',[1 0.5 0]);
% ylabel('blocklength');
% legend('opt \alpha_1','blocklength');
% 
% axis([-inf inf 0 2000]);
% 
% xline(dis_thred1,'-','Threshold1');
% xline(dis_thred2,'-.','Threshold2a');
% xline(dis_thred3,'--','Threshold2b');
% 
% 
% xlabel('Distance ratio');




figure (2)
syms D2
eqn_constraint1 = d1^-eta*eplsion1R - D2^-eta*eplsion2R == 0;
D2_sol1 = vpasolve(eqn_constraint1, D2);
eqn_constraint2 = (D2^-eta*eplsion2R*rho + 2)*d1^-eta*eplsion1R - (D2^-eta.*eplsion2R.*rho + 4)*D2^-eta*eplsion2R == 0;
D2_sol2 = vpasolve(eqn_constraint2, D2);

yyaxis left
plot(d2, real(opt_a1),'b');
hold on; grid on;
% plot(d2, real(opt_a1_mat(:,2)),'b','linewidth', 1.5);
% plot(d2, imag(opt_a1_mat(:,2)),'r','linewidth', 1.5);
yline(0.5,'-');

ylabel('\alpha_1');
axis([-inf inf 0 1]);

yyaxis right
plot(d2, real(opt_M), 'o', 'Color',[1 0.5 0]);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;
% plot(d2, real(opt_M_u1),'o','Color',[1 0.5 0]);
% plot(d2, real(opt_M_u2),'+','Color',[1 0.5 0]);
ylabel('blocklength');
% legend('opt \alpha_1','opt \alpha_1 by matlab real', 'opt \alpha_1 by matlab imag',...
%         'blocklength','blocklength by u1', 'blocklength by u2');

legend('opt \alpha_1','blocklength');
    
xline(double(D2_sol1(2)),'-','Threshold1');
xline(double(D2_sol2(2)),'-.','Threshold2');
% axis([-inf inf 0 2000]);


xlabel('user 2 distance');



