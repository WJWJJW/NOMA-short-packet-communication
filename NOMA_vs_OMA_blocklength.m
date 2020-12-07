clc; clear variables; close all;
N = 1e6;
N1 = 256;
N2 = 256;


eplsion1R = 10^-5;
eplsion2R = 10^-4;


BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

d1 = 400;
d2 = 900;
eta = 4;

h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

lamda1 = mean(abs(h1).^2);
lamda2 = mean(abs(h2).^2);


Pt = 35:1:40;                    %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

beta1 = 0.5;
beta2 = 1-beta1;

for u=1:length(Pt)
    rho = pt(u)/no;
    NOMA_opt_a1(u) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R))) ...
                     / (2*lamda2*eplsion2R*rho*(lamda1*eplsion1R-lamda2*eplsion2R));


    NOMA_opt_M(u) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R)))...
                     / (2*lamda2*eplsion2R)));
    
    % Equal power allocation
    rho1 = 0.5*pt(u) / (beta1*no);
    rho2 = 0.5*pt(u) / (beta2*no);
    OMA_opt_M1(u) = N1/(beta1*log2(1+(eplsion1R*lamda1*rho)));
    OMA_opt_M2(u) = N2/(beta2*log2(1+(eplsion2R*lamda2*rho)));
    OMA_opt_M(u) = OMA_opt_M1(u)+OMA_opt_M2(u);

end



figure(1)
plot(Pt, NOMA_opt_M, 'linewidth', 1.5);
hold on; grid on;
plot(Pt, OMA_opt_M, 'linewidth', 1.5);

xlabel('Transmit power (dBm)');
ylabel('M');
title('Required blocklength vs Transmit Power');
legend('NOMA','OMA');





