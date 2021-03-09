clc; clear variables; close all;
tic
%% Basic setting
N = 1e6;
eta = 4;
% User distance
d1 = 400;
d2 = 900;

% Rayleigh fading channel
h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

lambda1 = mean(abs(h1).^2);
lambda2 = mean(abs(h2).^2);

% AWGN noise
BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

% Parameter setting
M1 = 100;
M2 = 100;
N1 = 80;
N2 = 80;

w1 = 2^(N1/M1)-1;
w2 = 2^(N2/M2)-1;

Xi_1 = sqrt(1/(2*pi*(2^(2*N1/M1)-1)));
Xi_2 = sqrt(1/(2*pi*(2^(2*N2/M2)-1)));

eta_1 = w1 - (1/(2*sqrt(M1)*Xi_1));
tau_1 = w1 + (1/(2*sqrt(M1)*Xi_1));
eta_2 = w2 - (1/(2*sqrt(M2)*Xi_2));
tau_2 = w2 + (1/(2*sqrt(M2)*Xi_2));

% Transmit power in dBm
Pt = 40           
% Pt = -113:2:-74;
% Transmit power in linear scale
pt = (10^-3)*10^(Pt/10);

% Power allocaitn coefficient
alpha1 = 0.1:0.01:0.5;
alpha2 = 1 - alpha1;
%%

for aa=1:length(alpha1)
    rho = (10^-3)*db2pow(Pt)/ no;
    %% Q-function w/ NO approximation
    % Received SINR
    gamma2(aa,:) = (alpha2(aa)*pt.*abs(h2).^2)./(alpha1(aa)*pt.*abs(h2).^2+no);
    gamma12(aa,:) = (alpha2(aa)*pt.*abs(h1).^2)./(alpha1(aa)*pt.*abs(h1).^2+no);
    gamma11(aa,:) = (alpha1(aa)*pt.*abs(h1).^2)./no;
    
    % Q-function argument
    x_2a(aa,:) = log2 (1+gamma2(aa,:)).*M2 - N2;
    x_2b(aa,:) = log2 (exp(1))*sqrt((1-1./(1+gamma2(aa,:)).^2) .*M2);
    ratio2(aa,:) = x_2a(aa,:)./x_2b(aa,:);

    x_12a(aa,:) = log2 (1+gamma12(aa,:)).*M1 - N2;
    x_12b(aa,:) = log2 (exp(1))*sqrt((1-1./(1+gamma12(aa,:)).^2) .*M1);
    ratio12(aa,:) = x_12a(aa,:)./x_12b(aa,:);

    x_11a(aa,:) = log2 (1+gamma11(aa,:)).*M1 - N1;
    x_11b(aa,:) = log2 (exp(1))*sqrt((1-1./(1+gamma11(aa,:)).^2) .*M1);
    ratio11(aa,:) = x_11a(aa,:)./x_11b(aa,:);
    
    % Q-function
    q2(aa,:) = qfunc(ratio2(aa,:));
    q12(aa,:) = qfunc(ratio12(aa,:)); 
    q11(aa,:) = qfunc(ratio11(aa,:)); 
    
    %% BLER from approximation Q by MATLAB
    fun1 = @(x) 1-exp(x./(-lambda2*rho*(alpha2(aa)-alpha1(aa)*x)));
    fun2 = @(x) 1-exp(x./(-lambda1*rho*(alpha2(aa)-alpha1(aa)*x)));
    fun3 = @(x) 1-exp(x./(-lambda1*rho*alpha1(aa)));

    MAT_Epsilon2(aa) = Xi_2*sqrt(M2)*integral(fun1,eta_2,tau_2);
    MAT_Epsilon12(aa) = Xi_2*sqrt(M1)*integral(fun2,eta_2,tau_2);
    MAT_Epsilon11(aa) = Xi_1*sqrt(M1)*integral(fun3,eta_1,tau_1);
    MAT_Epsilon1(aa) = MAT_Epsilon12(aa) + MAT_Epsilon11(aa);
    
    %% BLER from approximation Q + Riemann integral
    epsilon2(aa) = 1 - exp(-w2/(lambda2*rho*(alpha2(aa)-alpha1(aa)*w2)));
    epsilon12(aa) = 1 - exp(-w2/(lambda1*rho*(alpha2(aa)-alpha1(aa)*w2)));
    epsilon11(aa) = 1 - exp(-w1/(lambda1*rho*alpha1(aa)));
    epsilon1(aa) = epsilon12(aa) + epsilon11(aa);
    %% BLER from approximation Q + Riemann integral + High SNR approximation
    epsilon2_high_SNR(aa) = w2/(lambda2*rho*(alpha2(aa)-alpha1(aa)*w2));
    epsilon12_high_SNR(aa) = w2/(lambda1*rho*(alpha2(aa)-alpha1(aa)*w1));
    epsilon11_high_SNR(aa) = w1/(lambda1*rho*alpha1(aa));
    epsilon1_high_SNR(aa) = epsilon12_high_SNR(aa) + epsilon11_high_SNR(aa);
    
end


figure (1)
semilogy(alpha1, mean(q2,2),'r', 'linewidth', 1.5);
hold on; grid on;
semilogy(alpha1, mean(q12,2)+(1-mean(q12,2)).*mean(q11,2),'b', 'linewidth', 1.5);
semilogy(alpha1, mean(q12,2),'g', 'linewidth', 1.5);
semilogy(alpha1, mean(q11,2),'k', 'linewidth', 1.5);

semilogy(alpha1, MAT_Epsilon2,'or', 'linewidth',1.5);
semilogy(alpha1, MAT_Epsilon1, 'ob','linewidth',1.5);
semilogy(alpha1, MAT_Epsilon12,'og', 'linewidth',1.5);
semilogy(alpha1, MAT_Epsilon11, 'ok','linewidth',1.5);

semilogy(alpha1, epsilon2,'+r');
semilogy(alpha1, epsilon1, '+b');
semilogy(alpha1, epsilon12,'+g');
semilogy(alpha1, epsilon11, '+k');

semilogy(alpha1, epsilon2_high_SNR,'*r');
semilogy(alpha1, epsilon1_high_SNR, '*b');
semilogy(alpha1, epsilon12_high_SNR,'*g');
semilogy(alpha1, epsilon11_high_SNR, '*k');

% axis([0 40 1e-6 1]);
title('BLER vs Power Allocation Coefficient');
xlabel('Power Allocation Coefficient');
ylabel('BLER');

legend('BLER2 eq4','BLER1 eq4','BLER12 eq4','BLER11 eq4',...
       'BLER2 Q-appro','BLER1 Q-appro','BLER12 Q-appro','BLER11 Q-appro',...
       'BLER2 Q+R-appro','BLER1 Q+R-appro','BLER12 Q+R-appro','BLER11 Q+R-appro',...
       'BLER2 Q+R+High SNR-appro','BLER1 Q+R+High SNR-appro','BLER12 Q+R+High SNR-appro','BLER11 Q+R+High SNR-appro');

figure (2)
semilogy(alpha1, mean(q2,2),'r', 'linewidth', 1.5);
hold on; grid on;
semilogy(alpha1, mean(q12,2)+(1-mean(q12,2)).*mean(q11,2),'b', 'linewidth', 1.5);


semilogy(alpha1, MAT_Epsilon2,'or', 'linewidth',1.5);
semilogy(alpha1, MAT_Epsilon1, 'ob','linewidth',1.5);


semilogy(alpha1, epsilon2,'+r');
semilogy(alpha1, epsilon1, '+b');


semilogy(alpha1, epsilon2_high_SNR,'*r');
semilogy(alpha1, epsilon1_high_SNR, '*b');


% axis([0 40 1e-6 1]);
title('BLER vs Power Allocation Coefficient');
xlabel('Power Allocation Coefficient');
ylabel('BLER');

legend('BLER2 eq4','BLER1 eq4',...
       'BLER2 Q-appro','BLER1 Q-appro',...
       'BLER2 Q+R-appro','BLER1 Q+R-appro',...
       'BLER2 Q+R+High SNR-appro','BLER1 Q+R+High SNR-appro');



toc