clc; clear variables; close all;
N = 1e6;
N1 = 80;
N2 = 80;


eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 10;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt/ no;

eta = 4;
dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);
dis_thred2 = (2*eplsion2R/eplsion1R)^(1/eta);


d2 = 300;
d1 = 10:10:290;

syms otp_a1

for dd=1:length(d1)
    
    h1 = sqrt(1/2*d1(dd)^-eta)*(randn(1,N)+1i*randn(1,N));
    h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

    lamda1 = mean(abs(h1).^2);
    lamda2 = mean(abs(h2).^2);
    
    xo2 = (1+eplsion2R*lamda2*rho)/(1+eplsion2R*lamda2*rho*otp_a1);
    xo1 = 1+otp_a1*rho*(lamda1*eplsion1R-lamda2*eplsion2R);

    eqn_opt = (xo2) == (xo1);
    opt_a1_mat(dd,:) = vpasolve(eqn_opt,otp_a1);
    opt_a1(dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R)))...
              / (2*lamda2*eplsion2R*rho*(lamda1*eplsion1R-lamda2*eplsion2R));
          
    opt_M_u1(dd) = N1/log2(1+opt_a1(dd)*rho*(lamda1*eplsion1R-lamda2*eplsion2R));
    opt_M_u2(dd) = N1/log2(1+opt_a1(dd)*rho*(lamda1*eplsion1R-lamda2*eplsion2R));
    opt_M(dd) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R)))...
                     / (2*lamda2*eplsion2R)));
end

figure (1)
ratio = d2./d1;
yyaxis left
plot(ratio, real(opt_a1),'ob');
hold on; grid on;
plot(ratio, real(opt_a1_mat(:,2)),'b','linewidth', 1.5);
plot(ratio, imag(opt_a1_mat(:,2)),'r','linewidth', 1.5);
yline(0.5,'-');

ylabel('\alpha_1');
axis([-inf inf 0 1]);

yyaxis right
plot(ratio, real(opt_M),'Color',[1 0.5 0]);
hold on; grid on;
plot(ratio, real(opt_M_u1),'o','Color',[1 0.5 0]);
plot(ratio, real(opt_M_u2),'+','Color',[1 0.5 0]);
ylabel('blocklength');
legend('opt \alpha_1','opt \alpha_1 by matlab real', 'opt \alpha_1 by matlab imag',...
        'blocklength','blocklength by u1', 'blocklength by u2');

axis([-inf inf 0 2000]);

xline(dis_thred1,'-','Threshold1');
xline(dis_thred2,'-.','Threshold2');



xlabel('Distance ratio');

% title('Distance Ratio vs alpha1');

figure (2)
yyaxis left
plot(d1, real(opt_a1),'ob');
hold on; grid on;
plot(d1, real(opt_a1_mat(:,2)),'b','linewidth', 1.5);
plot(d1, imag(opt_a1_mat(:,2)),'r','linewidth', 1.5);
yline(0.5,'-');

ylabel('\alpha_1');
axis([-inf inf 0 1]);

yyaxis right
plot(d1, real(opt_M),'Color',[1 0.5 0]);
hold on; grid on;
plot(d1, real(opt_M_u1),'o','Color',[1 0.5 0]);
plot(d1, real(opt_M_u2),'+','Color',[1 0.5 0]);
ylabel('blocklength');
legend('opt \alpha_1','opt \alpha_1 by matlab real', 'opt \alpha_1 by matlab imag',...
        'blocklength','blocklength by u1', 'blocklength by u2');

% axis([-inf inf 0 2000]);

xline(dis_thred1,'-','Threshold1');
xline(dis_thred2,'-.','Threshold2');



xlabel('Distance ratio');



