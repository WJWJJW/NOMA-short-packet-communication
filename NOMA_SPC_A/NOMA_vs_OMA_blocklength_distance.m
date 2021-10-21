% This script analysis the relation between blocklength , power allocation
% coefficient and user 2 distance
% calculate blocklength and power allocation coefficient according
% different conditon

clc; clear variables; close all;
N = 1e6;
N1 = 256;
N2 = 256;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 30;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

% AWGN noise
% BW = 10^6;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt/ no;

eta = 4;

d1 = 100;

d2 = 101:1:300;
check1 = zeros(length(d2),1);
check2 = zeros(length(d2),1);


beta1 = 1/2;
beta2 = 1-beta1;

% 
% beta1 = 1;
% beta2 = 1;

h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
lamda1 = mean(abs(h1).^2);
[delta] = delta_finder (lamda1*rho*eplsion1R);

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
    % Blocklength for OMA users
    % Equal power allocation
    rho1 = (0.5*pt) / (beta1*no);
    rho2 = (0.5*pt) / (beta2*no);
    OMA_opt_M1(dd) = N1/(beta1*log2(1+(eplsion1R*lamda1*rho1)));
    OMA_opt_M2(dd) = N2/(beta2*log2(1+(eplsion2R*lamda2*rho2)));
        
    OMA_opt_M(dd) = OMA_opt_M1(dd)+OMA_opt_M2(dd);


end

% min_lamda2 = (lamda1^2 * 4 * (eplsion1R+eplsion2R)^2 / (delta*eplsion1R*eplsion2R^2*rho))^(1/3);
% f1 = sqrt(1+0.5*rho*lamda1*eplsion1R);
% ff = (-1+f1)/f1;
% min_lamda2 = (2*lamda1*eplsion1R*ff) / (eplsion2R*(lamda1*eplsion1R*rho - 2*ff));
% min_d2 = min_lamda2 ^(-1/eta)
% 
% g = (-1+sqrt(1+4*delta^2*(1-delta)*lamda1*eplsion1R*rho))/(2*delta);
% min_lamda22 = (lamda1*eplsion1R*g)/((0.5*lamda1*eplsion1R*rho-g)*eplsion2R);
% min_d22 = min_lamda22 ^(-1/eta)

figure (1)
syms D2;
eqn_constraint1 = d1^-eta*eplsion1R - D2^-eta*eplsion2R == 0;
D2_sol1 = vpasolve(eqn_constraint1, D2);
eqn_constraint2 = (D2^-eta*eplsion2R*rho + 2)*d1^-eta*eplsion1R - (D2^-eta.*eplsion2R.*rho + 4)*D2^-eta*eplsion2R == 0;
D2_sol2 = vpasolve(eqn_constraint2, D2);

% eqn_test1 = (beta1*eplsion1R*lamda1*rho1+beta2*eplsion2R*(D2)^(-eta)*rho2)*g ...
%     - (beta1*(eplsion1R*lamda1*rho1)/(1+eplsion1R*lamda1*rho1))*(beta2*(eplsion2R*(D2)^(-eta)*rho2)/(1+eplsion2R*(D2)^(-eta)*rho2));

% D2_sol222 = vpasolve(eqn_test1, D2)

% eqn_OMA = 1/(beta1*log2(1+(eplsion1R*lamda1*rho1))) + 1/(beta2*log2(1+(eplsion2R*(D2)^(-eta)*rho2))) ...
%         - 1/log2(1+((-1+sqrt(1+4*delta^2*rho*lamda1*eplsion1R*(1-delta))) / (2*delta))) == 0;


% eqn_OMA = (beta1*log(1+(eplsion1R*lamda1*rho1))+beta2*log(1+(eplsion2R*(D2)^(-eta)*rho2)))...
%         *(log(1+((-1+sqrt(1+4*delta^2*rho*lamda1*eplsion1R*(1-delta))) / (2*delta))))...
%         - beta1*log(1+(eplsion1R*lamda1*rho1))*beta2*log(1+(eplsion2R*(D2)^(-eta)*rho2)) == 0


% eqn_OMA = (eplsion1R*lamda1 + eplsion2R*(D2)^(-eta)) * ((-1+sqrt(1+4*delta^2*rho*lamda1*eplsion1R*(1-delta))) / (2*delta)) ...
%     -0.5*rho*eplsion1R*lamda1*eplsion2R*(D2)^(-eta) == 0;

% eqn_OMA = 4*delta*(1-delta)*(eplsion1R*lamda1 + eplsion2R*(D2)^(-eta))^2 == ...
%             eplsion2R*(D2)^(-eta)*(delta*rho*eplsion1R*lamda1*eplsion2R*(D2)^(-eta)+2*(eplsion1R*lamda1 + eplsion2R*(D2)^(-eta)));

% D2_sol3 = abs(vpasolve(eqn_OMA, D2))

term1 = beta1*log(1+eplsion1R*lamda1*rho1);
term2 = log(1+((-1+sqrt(1+4*delta^2*rho*lamda1*eplsion1R*(1-delta))) / (2*delta)))
d2_min_thred = ((exp((term1*term2)/((term1-term2)*beta2))-1) / (eplsion2R*rho2))^(-1/eta)


plot(d2, opt_M, 'b','linewidth',1.5);
hold on; grid on;
plot(d2, OMA_opt_M, 'r','linewidth',1.5);

ylabel('Blocklength (Channel uses)');


legend('NOMA: Blocklength by (6.7) or (6.12)','OMA: Blocklength by Sum of (6.14)');
    
% xline(double(D2_sol1(2)),'-','Threshold1');
xx = xline(double(D2_sol2(2)),'-.','Bound (6.9)','HandleVisibility','off');
xx.FontName = 'Times New Roman';

% xline(double(min_d2),'-.','Bound1');
% xline(double(min_d22),'-.','Bound2');


% xline(double(D2_sol1(2)),'-','D2 asym');
xxx = xline(d2_min_thred,'-','Performance Bound by Thm. 4','HandleVisibility','off');
xxx.FontName = 'Times New Roman';
xlabel('Distance of far user (meter)');

set(gca, 'FontName', 'Times New Roman');

figure(2)

plot(d2, opt_M, 'b','linewidth',1.5);
hold on; grid on;
plot(d2, OMA_opt_M, 'r','linewidth',1.5);
plot(d2, OMA_opt_M1, 'g','linewidth',1.5);
plot(d2, OMA_opt_M2, 'm','linewidth',1.5);

xx = xline(double(D2_sol2(2)),'-.','Bound (6.9)','HandleVisibility','off');
xx.FontName = 'Times New Roman';

% xline(double(min_d2),'-.','Bound1');
% xline(double(min_d22),'-.','Bound2');


% xline(double(D2_sol1(2)),'-','D2 asym');
xxx = xline(d2_min_thred,'-','Performance Bound by Thm. 4','HandleVisibility','off');
xxx.FontName = 'Times New Roman';

legend('NOMA: Blocklength by (6.7) or (6.12)','OMA: Blocklength by Sum of (6.14)','OMA near user', 'OMA far user');

xlabel('Distance of far user (meter)');
ylabel('Blocklength (Channel uses)');

set(gca, 'FontName', 'Times New Roman');