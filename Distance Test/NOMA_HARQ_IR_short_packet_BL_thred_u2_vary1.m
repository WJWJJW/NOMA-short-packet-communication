% This script analysis the relation between blocklength , power allocation
% coefficient and user 2 distance
% calculate blocklength and power allocation coefficient according
% different conditon

clc; clear variables; close all;
N = 1e7;
N1 = 256;
N2 = 256;

Eplsion1R = [10^-5 10^-4];
Eplsion2R = [10^-4 10^-5];

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
d2 = 100:2:300;


h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
lamda1 = mean(abs(h1).^2);

opt_a1 = zeros(length(d2),length(Eplsion1R),length(Eplsion2R));
opt_M = zeros(length(d2),length(Eplsion1R),length(Eplsion2R));

for er1 = 1:length(Eplsion1R)
    eplsion1R = Eplsion1R(er1);
    [delta] = delta_finder (lamda1*rho*eplsion1R);
    for er2 = 1:length(Eplsion2R)
        eplsion2R = Eplsion2R(er2);
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
            opt_a1(dd,er1,er2) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
                      4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m))) ...
                      / (2*lamda2*eplsion2R_m*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m));
            opt_M(dd,er1,er2) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)))...
                                 / (2*lamda2*eplsion2R_m)));
        end
    end
end

D2_sol2 = zeros(length(Eplsion1R),length(Eplsion2R));
for er1 = 1:length(Eplsion1R)
    for er2 = 1:length(Eplsion2R)
        eplsion1R = Eplsion1R(er1);
        eplsion2R = Eplsion2R(er2);
        syms D2
        eqn_constraint2 = (D2^-eta*eplsion2R*rho + 2)*d1^-eta*eplsion1R - (D2^-eta.*eplsion2R.*rho + 4)*D2^-eta*eplsion2R == 0;
        d2_sol2 = vpasolve(eqn_constraint2, D2);
        D2_sol2(er1,er2) = d2_sol2(2);
    end
end


figure (1)
subplot(2,2,1);
yyaxis left
plot(d2, real(opt_a1(:,1,1)),'b', 'linewidth',1.5);
hold on; grid on;

yline(0.5,'-','HandleVisibility','off');

ylabel('\alpha_1');
% axis([-inf inf 0 0.5]);

yyaxis right
plot(d2, real(opt_M(:,1,1)), 'Color',[1 0.5 0], 'linewidth',1.5);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('Blocklength (Channel uses)');
legend('\alpha_1^{opt}','M_{opt}^{1 pair}');
    
% xx = xline(double(D2_sol2(2)),'-.','Bound (38)','HandleVisibility','off');
% xx.FontName = 'Times New Roman';

xlabel('Distance of far user (meter)');
set(gca, 'FontName', 'Times New Roman');

subplot(2,2,2);
yyaxis left
plot(d2, real(opt_a1(:,1,2)),'b', 'linewidth',1.5);
hold on; grid on;

yline(0.5,'-','HandleVisibility','off');

ylabel('\alpha_1');
% axis([-inf inf 0 0.5]);

yyaxis right
plot(d2, real(opt_M(:,1,2)), 'Color',[1 0.5 0], 'linewidth',1.5);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('Blocklength (Channel uses)');
legend('\alpha_1^{opt}','M_{opt}^{1 pair}');
    
% xx = xline(double(D2_sol2(2)),'-.','Bound (38)','HandleVisibility','off');
% xx.FontName = 'Times New Roman';

xlabel('Distance of far user (meter)');
set(gca, 'FontName', 'Times New Roman');

subplot(2,2,3);
yyaxis left
plot(d2, real(opt_a1(:,2,1)),'b', 'linewidth',1.5);
hold on; grid on;

yline(0.5,'-','HandleVisibility','off');

ylabel('\alpha_1');
% axis([-inf inf 0 0.5]);

yyaxis right
plot(d2, real(opt_M(:,2,1)), 'Color',[1 0.5 0], 'linewidth',1.5);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('Blocklength (Channel uses)');
legend('\alpha_1^{opt}','M_{opt}^{1 pair}');
    
% xx = xline(double(D2_sol2(2)),'-.','Bound (38)','HandleVisibility','off');
% xx.FontName = 'Times New Roman';

xlabel('Distance of far user (meter)');
set(gca, 'FontName', 'Times New Roman');
subplot(2,2,4);
yyaxis left
plot(d2, real(opt_a1(:,2,2)),'b', 'linewidth',1.5);
hold on; grid on;

yline(0.5,'-','HandleVisibility','off');

ylabel('\alpha_1');
% axis([-inf inf 0 0.5]);

yyaxis right
plot(d2, real(opt_M(:,2,2)), 'Color',[1 0.5 0], 'linewidth',1.5);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('Blocklength (Channel uses)');
legend('\alpha_1^{opt}','M_{opt}^{1 pair}');
    
% xx = xline(double(D2_sol2(2)),'-.','Bound (38)','HandleVisibility','off');
% xx.FontName = 'Times New Roman';

xlabel('Distance of far user (meter)');
set(gca, 'FontName', 'Times New Roman');




figure (2)
yyaxis left
plot(d2, real(opt_a1(:,1,1)),'b', 'linewidth',1.5);
hold on; grid on;

yline(0.5,'-','HandleVisibility','off');

ylabel('\alpha_1');
% axis([-inf inf 0 0.5]);

yyaxis right
plot(d2, real(opt_M(:,1,1)), 'Color',[1 0.5 0], 'linewidth',1.5);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('Blocklength (Channel uses)');
legend('\alpha_1^{opt}','M_{opt}^{1 pair}','FontSize',14);
    
xx = xline(double(D2_sol2(1,1)),'-.','Bound (38)','HandleVisibility','off');
xx.FontName = 'Times New Roman';
xx.FontSize = 14;

xlabel('Distance of far user (meter)');
set(gca, 'FontName', 'Times New Roman','FontSize',16);

figure (3)
yyaxis left
plot(d2, real(opt_a1(:,1,2)),'b', 'linewidth',1.5);
hold on; grid on;

yline(0.5,'-','HandleVisibility','off');

ylabel('\alpha_1');
% axis([-inf inf 0 0.5]);

yyaxis right
plot(d2, real(opt_M(:,1,2)), 'Color',[1 0.5 0], 'linewidth',1.5);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('Blocklength (Channel uses)');
legend('\alpha_1^{opt}','M_{opt}^{1 pair}','FontSize',14);
    
xx = xline(double(D2_sol2(1,2)),'-.','Bound (38)','HandleVisibility','off');
xx.FontName = 'Times New Roman';
xx.FontSize = 14;

xlabel('Distance of far user (meter)');
set(gca, 'FontName', 'Times New Roman','FontSize',16);

figure (4)
yyaxis left
plot(d2, real(opt_a1(:,2,1)),'b', 'linewidth',1.5);
hold on; grid on;

yline(0.5,'-','HandleVisibility','off');

ylabel('\alpha_1');
% axis([-inf inf 0 0.5]);

yyaxis right
plot(d2, real(opt_M(:,2,1)), 'Color',[1 0.5 0], 'linewidth',1.5);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('Blocklength (Channel uses)');
legend('\alpha_1^{opt}','M_{opt}^{1 pair}','FontSize',14);
    
xx = xline(double(D2_sol2(2,1)),'-.','Bound (38)','HandleVisibility','off');
xx.FontName = 'Times New Roman';
xx.FontSize = 14;

xlabel('Distance of far user (meter)');
set(gca, 'FontName', 'Times New Roman','FontSize',16);

figure (5)
yyaxis left
plot(d2, real(opt_a1(:,2,2)),'b', 'linewidth',1.5);
hold on; grid on;

yline(0.5,'-','HandleVisibility','off');

ylabel('\alpha_1');
% axis([-inf inf 0 0.5]);

yyaxis right
plot(d2, real(opt_M(:,2,2)), 'Color',[1 0.5 0], 'linewidth',1.5);
% set(gca, 'YTick', 120:1:125, 'YTickLabel', 120:1:125);
hold on; grid on;

ylabel('Blocklength (Channel uses)');
legend('\alpha_1^{opt}','M_{opt}^{1 pair}','FontSize',14);
    
xx = xline(double(D2_sol2(2,2)),'-.','Bound (38)','HandleVisibility','off');
xx.FontName = 'Times New Roman';
xx.FontSize = 14;

xlabel('Distance of far user (meter)');
set(gca, 'FontName', 'Times New Roman','FontSize',16);


