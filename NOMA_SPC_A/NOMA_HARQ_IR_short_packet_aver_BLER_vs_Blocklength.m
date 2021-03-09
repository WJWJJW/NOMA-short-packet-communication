% This script analysis relation between blocklength and BLER
% Find optimal blocklength first
% Show BLER against blocklength

clc; clear variables; close all;
N = 1e6;
N1 = 80;
N2 = 80;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 20;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)


rho = pt/ no;


d1 = 100;
d2 = 150;
eta = 4;

h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

lamda1 = mean(abs(h1).^2);
lamda2 = mean(abs(h2).^2);

contraint = (lamda1*eplsion1R)/(lamda2*eplsion2R);
[delta] = delta_finder (lamda1*rho*eplsion1R);
theta = contraint*delta;

check = (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - (lamda2*eplsion2R*rho + 4)*lamda2*eplsion2R > 0;
switch (check)
    case 0
        eplsion2R_m = theta*eplsion2R;
    case 1
        eplsion2R_m = eplsion2R;
end
oopt_M = N1/log2(1+((-lamda1*eplsion1R + ...
                sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)))...
                / (2*lamda2*eplsion2R_m)));

opt_a1 = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
              4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m))) ...
              / (2*lamda2*eplsion2R_m*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m));
a = opt_a1;
MM = 100:10:4000;

for m = 1:length(MM)
    BLER12(m) = (2^(N2/ MM(m))-1)/(lamda1*rho*(1-a-a*(2^(N2/ MM(m))-1)));
    BLER11(m) = (2^(N1/MM(m))-1)/(lamda1*a*rho);
    BLER1(m) = (2^(N2/ MM(m))-1)/(lamda1*rho*(1-a-a*(2^(N2/ MM(m))-1)))+(2^(N1/MM(m))-1)/(lamda1*a*rho);
    BLER2(m) = (2^(N2/MM(m))-1)/(lamda2*rho*(1-a-a*(2^(N2/MM(m))-1)));
    
end


figure (1)
semilogy(MM, BLER2, 'r', 'linewidth', 1.5);
hold on; grid on;
semilogy(MM, BLER1, 'b', 'linewidth', 1.5);
semilogy(MM, BLER11, 'm', 'linewidth', 1.5);
semilogy(MM, BLER12, 'g', 'linewidth', 1.5);

xlabel('Blocklength');
ylabel('BLER');
title('BLER vs Blocklength');
legend('Far user','Near user','Near user11','Near user12');
xline(oopt_M,'-');
yline(1e-4,'-');
yline(1e-5,'-');



