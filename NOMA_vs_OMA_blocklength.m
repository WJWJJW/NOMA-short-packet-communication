% This script analysis performance between NOMA OMA
% 2 user case
clc; clear variables; close all;
N = 1e6;
N1 = 80;
N2 = 80;


eplsion1R = 10^-5;
eplsion2R = 10^-4;


BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

d1 = 100;
d2 = 150;
eta = 4;

h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

% h1 = sqrt(1/2*(1+d1^eta)^-1)*(randn(1,N)+1i*randn(1,N));
% h2 = sqrt(1/2*(1+d2^eta)^-1)*(randn(1,N)+1i*randn(1,N));


lamda1 = mean(abs(h1).^2);
lamda2 = mean(abs(h2).^2);


Pt = 0:1:20;                    %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

beta1 = 0.5;
beta2 = 1-beta1;

delta = 0.5;

NOMA_opt_M = zeros(1,length(Pt));
OMA_opt_M1 = zeros(1,length(Pt));
OMA_opt_M2 = zeros(1,length(Pt));
OMA_opt_M = zeros(1,length(Pt));


for u=1:length(Pt)
    Rho(u) = Pt(u)-No;
    
    rho = pt(u)/no;
    rho*lamda1*eplsion1R;
    
    % Blocklength for NOMA users
    contraint = (lamda1*eplsion1R)/(lamda2*eplsion2R);
    theta = contraint*delta;
        
    check = (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - (lamda2*eplsion2R.*rho + 4)*lamda2*eplsion2R > 0
    switch (check)
        case 0
            eplsion2R_m = theta*eplsion2R;
        case 1
            eplsion2R_m = eplsion2R;
    end
    opt_a1(u) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
                  4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m))) ...
                  / (2*lamda2*eplsion2R_m*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m));
              
    opt_a11(u) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
                  4*(lamda2*eplsion2R_m)^2)) ...
                  / (2*lamda2*eplsion2R_m);
    
    opt_a1(u)*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)
              
    NOMA_opt_M(u) = N1/log2(1+((-lamda1*eplsion1R + ...
                    sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)))...
                    / (2*lamda2*eplsion2R_m)));
    
                
    NOMA_opt_M1(u) = N1/(opt_a1(u)*rho*(lamda1*eplsion1R - lamda2*eplsion2R_m));
    NOMA_opt_M21(u) = N1/log2((1+rho*lamda2*eplsion2R_m)/(1+rho*lamda2*eplsion2R_m*opt_a1(u)));
    NOMA_opt_M22(u) = N1/log2(1+opt_a11(u)*rho*(lamda1*eplsion1R - lamda2*eplsion2R_m));
    

    % Blocklength for OMA users
    % Equal power allocation
    rho1 = (0.5*pt(u)) / (beta1*no);
    rho2 = (0.5*pt(u)) / (beta2*no);
    OMA_opt_M1(u) = N1/(beta1*log2(1+(eplsion1R*lamda1*rho1)));
    OMA_opt_M2(u) = N2/(beta2*log2(1+(eplsion2R*lamda2*rho2)));
    
    OMA_opt_M1_appro(u) = N1/(beta1*(eplsion1R*lamda1*rho1));
    OMA_opt_M2_appro(u) = N2/(beta2*(eplsion2R*lamda2*rho2));
    
    OMA_opt_M(u) = OMA_opt_M1(u)+OMA_opt_M2(u);
    OMA_opt_M_appro(u) = OMA_opt_M1_appro(u)+OMA_opt_M2_appro(u);

end



figure(1)
plot(Pt, NOMA_opt_M, 'b', 'linewidth', 1.5);
hold on; grid on;
plot(Pt, OMA_opt_M, 'r', 'linewidth', 1.5);
plot(Pt, OMA_opt_M_appro, 'g', 'linewidth', 1.5);

xlabel('Transmit power (dBm)');
ylabel('M');
title('Required blocklength vs Transmit Power');
legend('NOMA','OMA');



figure(2)
plot(Pt, NOMA_opt_M, 'b', 'linewidth', 1.5);
hold on; grid on;
plot(Pt, OMA_opt_M1, '-or', 'linewidth', 1.5);
plot(Pt, OMA_opt_M2, '-sr', 'linewidth', 1.5);
plot(Pt, OMA_opt_M - NOMA_opt_M, 'g', 'linewidth', 1.5);

xlabel('Transmit power (dBm)');
ylabel('M');
title('Required blocklength vs Transmit Power');
legend('NOMA','OMA1','OMA2','Gap');



figure(3)
plot(Pt, NOMA_opt_M, 'b', 'linewidth', 1.5);
hold on; grid on;
plot(Pt, NOMA_opt_M21, '*b', 'linewidth', 1.5);
plot(Pt, NOMA_opt_M22, 'sb', 'linewidth', 1.5);

plot(Pt, NOMA_opt_M1, 'ob', 'linewidth', 1.5);

plot(Pt, OMA_opt_M1, '-or', 'linewidth', 1.5);
plot(Pt, OMA_opt_M2, '-sr', 'linewidth', 1.5);

plot(Pt, OMA_opt_M1_appro, '-og', 'linewidth', 1.5);
plot(Pt, OMA_opt_M2_appro, '-sg', 'linewidth', 1.5);





