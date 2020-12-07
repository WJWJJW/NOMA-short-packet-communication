clc; clear variables; close all;

N = 40;
M = 100;
R = N/M;
Rho=0:0.01:40;
rho=10^-3.*db2pow(Rho);

Xi = sqrt(1/(2*pi*(2^2*R-1)));
Omega = 2^R - 1;
Nu = Omega - 1/(2*sqrt(M)*Xi);
Tau = Omega + 1/(2*sqrt(M)*Xi);

for u = 1:length(rho)
    Q(u) = qfunc((M*log2(1+rho(u))-N) / (log2(exp(1))*sqrt(M*(1-1/(1+rho(u))^2))));

    if rho(u) <= Nu
        Q_appro(u) = 1;
    elseif rho(u) > Nu && rho(u) < Tau
        Q_appro(u) = 0.5 - Xi*sqrt(M)*(rho(u)-Omega);
    else
        Q_appro(u) = 0;
    end
end

figure(1)
plot(Rho, Q, 'b', 'linewidth', 1.5);
hold on;grid on;
plot(Rho, Q_appro, 'r', 'linewidth', 1.5);

xlabel('SNR(dB)');
ylabel('Q(.)');
title('Illustration of Q-function linearization');
legend('Q-function','Linearization');

