clc; clear variables; close all;
N = 1e6;
N1 = 80;
N2 = 80;


eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 30;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)


% rho = pt/ no
rho = db2pow(90);


d1 = 10;
d2 = 15;
eta = 4;

h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

lamda1 = mean(abs(h1).^2);
lamda2 = mean(abs(h2).^2);

t_delta1 = (lamda1*eplsion1R*(lamda2*eplsion2R*rho+2))/(2*lamda2*eplsion2R)-1
t_delta2 = t_delta1/100


if d2/d1 > (eplsion2R/eplsion1R)^0.25
oopt_a1 = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
          4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R))) ...
          / (2*lamda2*eplsion2R*rho*(lamda1*eplsion1R-lamda2*eplsion2R));
oopt_M = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R)))...
                     / (2*lamda2*eplsion2R)));
else
oopt_a1 = (-lamda1*eplsion1R - (t_delta2+1)*lamda2*eplsion2R + ...
          sqrt((lamda1*eplsion1R + (t_delta2+1)*lamda2*eplsion2R)^2 ...
          + 4*lamda1*eplsion1R*(lamda2*eplsion2R)^2*(t_delta2+1)*rho)) ...
          / (2*lamda1*eplsion1R*lamda2*eplsion2R*rho)
oopt_M = N1/log2(1+(lamda1*eplsion1R*oopt_a1*rho)/(t_delta2+1))
end

a = 0.3;
MM = 100:10:4000;

for m = 1:length(MM)
    BLER12(m) = (2^(N2/ MM(m))-1)/(lamda1*rho*(1-a-a*(2^(N2/ MM(m))-1)));
    BLER11(m) = (2^(N1/MM(m))-1)/(lamda1*a*rho);
    BLER1(m) = (2^(N2/ MM(m))-1)/(lamda1*rho*(1-a-a*(2^(N2/ MM(m))-1)))+(2^(N1/MM(m))-1)/(lamda1*a*rho);
    BLER2(m) = (2^(N2/MM(m))-1)/(lamda2*rho*(1-a-a*(2^(N2/MM(m))-1)));
    
    BLER1A(m) = (t_delta2+1)*(2^(N1/MM(m))-1)/(lamda1*a*rho);
end


figure (1)
semilogy(MM, BLER2, 'r', 'linewidth', 1.5);
hold on; grid on;
semilogy(MM, BLER1, 'b', 'linewidth', 1.5);
semilogy(MM, BLER11, 'm', 'linewidth', 1.5);
semilogy(MM, BLER12, 'g', 'linewidth', 1.5);
semilogy(MM, BLER1A, 'ok');
xlabel('Blocklength');
ylabel('BLER');
title('BLER vs Blocklength');
legend('Far user','Near user','Near user11','Near user12');
xline(oopt_M,'-');
yline(1e-4,'-');
yline(1e-5,'-');



