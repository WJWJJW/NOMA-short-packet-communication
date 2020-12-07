clc; clear variables; close all;
tic
%% Basic setting
N = 1e6;
eta = 4;
% User distance
d1 = 400;
d2 = 1700;

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
M1 = 300;
M2 = 300;
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
Pt = 0:2:40;                
% Pt = -113:2:-74;
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);

% Power allocaitn coefficient
alpha1 = 0.1;
alpha2 = 1 - alpha1;
%%

for u=1:length(Pt)
    rho = (10^-3)*db2pow(Pt(u))/ no;
    %% Q-function w/ NO approximation
    % Received SINR
    gamma2(u,:) = (alpha2*pt(u).*abs(h2).^2)./(alpha1*pt(u).*abs(h2).^2+no);
    gamma12(u,:) = (alpha2*pt(u).*abs(h1).^2)./(alpha1*pt(u).*abs(h1).^2+no);
    gamma11(u,:) = (alpha1*pt(u).*abs(h1).^2)./no;
    
    % Q-function argument
    x_2a(u,:) = log2 (1+gamma2(u,:)).*M2 - N2;
    x_2b(u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma2(u,:)).^2) .*M2);
    ratio2(u,:) = x_2a(u,:)./x_2b(u,:);

    x_12a(u,:) = log2 (1+gamma12(u,:)).*M2 - N2;
    x_12b(u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma12(u,:)).^2) .*M2);
    ratio12(u,:) = x_12a(u,:)./x_12b(u,:);

    x_11a(u,:) = log2 (1+gamma11(u,:)).*M1 - N1;
    x_11b(u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma11(u,:)).^2) .*M1);
    ratio11(u,:) = x_11a(u,:)./x_11b(u,:);
    
    % Q-function
    q2(u,:) = qfunc(ratio2(u,:));
    q12(u,:) = qfunc(ratio12(u,:)); 
    q11(u,:) = qfunc(ratio11(u,:)); 
    
    %% BLER from approximation Q by MATLAB
    fun1 = @(x) 1-exp(x./(-lambda2*rho*(alpha2-alpha1*x)));
    fun2 = @(x) 1-exp(x./(-lambda1*rho*(alpha2-alpha1*x)));
    fun3 = @(x) 1-exp(x./(-lambda1*rho*alpha1));

    MAT_Epsilon2(u) = Xi_2*sqrt(M2)*integral(fun1,eta_2,tau_2);
    MAT_Epsilon12(u) = Xi_2*sqrt(M2)*integral(fun2,eta_2,tau_2);
    MAT_Epsilon11(u) = Xi_1*sqrt(M1)*integral(fun3,eta_1,tau_1);
    MAT_Epsilon1(u) = MAT_Epsilon12(u) + MAT_Epsilon11(u);
    
    %% BLER from approximation Q + Riemann integral
    epsilon2(u) = 1 - exp(-w2/(lambda2*rho*(alpha2-alpha1*w2)));
    epsilon12(u) = 1 - exp(-w2/(lambda1*rho*(alpha2-alpha1*w2)));
    epsilon11(u) = 1 - exp(-w1/(lambda1*rho*alpha1));
    epsilon1(u) = epsilon12(u) + epsilon11(u);
    %% BLER from approximation Q + Riemann integral + High SNR approximation
    epsilon2_high_SNR(u) = w2/(lambda2*rho*(alpha2-alpha1*w2));
    epsilon12_high_SNR(u) = w2/(lambda1*rho*(alpha2-alpha1*w2));
    epsilon11_high_SNR(u) = w1/(lambda1*rho*alpha1);
    epsilon1_high_SNR(u) = epsilon12_high_SNR(u) + epsilon11_high_SNR(u);
    
end

%% Plot
% figure (1)
% plot(Pt, mean(q2,2),'r', 'linewidth', 1.5);
% hold on; grid on;
% plot(Pt, mean(q12,2)+(1-mean(q12,2)).*mean(q11,2),'b', 'linewidth', 1.5);
% plot(Pt, mean(q12,2),'og', 'linewidth', 1.5);
% plot(Pt, mean(q11,2),'ok', 'linewidth', 1.5);
% 
% title('BLER vs Transmit SNR');
% xlabel('Transmit SNR (Rho in dBm)');
% ylabel('BLER');
% 
% legend('epsilon2','epsilon1','epsilon12','epsilon11');

figure (1)
semilogy(Pt, mean(q2,2),'r', 'linewidth', 1.5);
hold on; grid on;
semilogy(Pt, mean(q12,2)+(1-mean(q12,2)).*mean(q11,2),'b', 'linewidth', 1.5);
semilogy(Pt, mean(q12,2),'g', 'linewidth', 1.5);
semilogy(Pt, mean(q11,2),'k', 'linewidth', 1.5);

semilogy(Pt, MAT_Epsilon2,'or', 'linewidth',1.5);
semilogy(Pt, MAT_Epsilon1, 'ob','linewidth',1.5);
semilogy(Pt, MAT_Epsilon12,'og', 'linewidth',1.5);
semilogy(Pt, MAT_Epsilon11, 'ok','linewidth',1.5);

semilogy(Pt, epsilon2,'+r');
semilogy(Pt, epsilon1, '+b');
semilogy(Pt, epsilon12,'+g');
semilogy(Pt, epsilon11, '+k');

semilogy(Pt, epsilon2_high_SNR,'*r');
semilogy(Pt, epsilon1_high_SNR, '*b');
semilogy(Pt, epsilon12_high_SNR,'*g');
semilogy(Pt, epsilon11_high_SNR, '*k');

% axis([0 40 1e-6 1]);
title('BLER vs Transmitted power');
xlabel('Transmitted power (in dBm)');
ylabel('BLER');
dim = [0.15 0.5 0.005 0.01];
str = {'Line : No approximation',...
       'Circle : Q-function linearization',...
       'Cross : Q-function linearization + ','Riemann integral',...
       'Star : Q-function linearization + ','Riemann integral + ','High SNR approximation'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
legend('BLER2','BLER1','BLER12','BLER11');
   
figure (2)
semilogy(Pt, mean(q2,2),'r', 'linewidth', 1.5);
hold on; grid on;
semilogy(Pt, mean(q12,2)+(1-mean(q12,2)).*mean(q11,2),'b', 'linewidth', 1.5);


semilogy(Pt, MAT_Epsilon2,'or', 'linewidth',1.5);
semilogy(Pt, MAT_Epsilon1, 'ob','linewidth',1.5);


semilogy(Pt, epsilon2,'+r');
semilogy(Pt, epsilon1, '+b');


semilogy(Pt, epsilon2_high_SNR,'*r');
semilogy(Pt, epsilon1_high_SNR, '*b');


title('BLER vs Transmitted power');
xlabel('Transmitted power (in dBm)');
ylabel('BLER');
dim = [0.15 0.5 0.005 0.01];
str = {'Line : No approximation',...
       'Circle : Q-function linearization',...
       'Cross : Q-function linearization + ','Riemann integral',...
       'Star : Q-function linearization + ','Riemann integral + ','High SNR approximation'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

legend('BLER2','BLER1');


toc