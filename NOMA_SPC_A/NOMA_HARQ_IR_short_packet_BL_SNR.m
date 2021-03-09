clc; clear variables; close all;
N1=80;
N2=80;
eplsion1R=10^-5;
eplsion2R=10^-4;
% rho = 10^(RHO_db/10);
% rho = (10^-3)*10^(RHO_db/10);

d1 = 400;
d2 = 900;
eta = 4;


h1 = sqrt(1/2*d1^-eta)*(randn(1,10000)+1i*randn(1,10000));
h2 = sqrt(1/2*d2^-eta)*(randn(1,10000)+1i*randn(1,10000));

lamda1 = mean(abs(h1).^2);
lamda2 = mean(abs(h2).^2);

% lamda1 = 1/(1+d1^eta);
% lamda2 = 1/(1+d2^eta);

BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)


Pt = 50:55;
pt = (10^-3).* db2pow(Pt);
rho = pt ./ no;

alpha1 = (eplsion2R * lamda2) / (eplsion1R * lamda1 - eplsion2R * lamda2)
for u = 1:length(rho)
    x2 = (1+eplsion2R*lamda2*rho(u))/(1+eplsion2R*lamda2*rho(u)*alpha1);
    M2(u) = N2/(log2 (x2));

    x1 = 1+alpha1*rho(u)*(lamda1*eplsion1R-lamda2*eplsion2R);
    M1(u) = N1/(log2 (x1));

end



figure;
plot(Pt, M1, 'ob', 'linewidth', 1.5);
hold on; grid on;
plot(Pt, M2, 'r', 'linewidth', 1.5);


xlabel('SNR');
ylabel('M');
title('Required blocklength vs Transmitted SNR');
legend('$M_1$','$M_2$', 'Interpreter', 'latex');



