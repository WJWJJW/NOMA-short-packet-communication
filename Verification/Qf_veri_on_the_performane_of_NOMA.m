clear all;close all;clc;

tic

N = 1e6

% Transmit power in dB
RHO = 0:2:50;                
% Transmit power in linear scale
rho = 10.^(RHO/10);

% Power allocaitn coefficient
alpha1 = [0.1 0.2];
alpha2 = 1 - alpha1;

% Rayleigh fading channel
lambda = 1;


% Instantaneous BLER by Q-function
M1 = 100;
M2 = 100;
N1 = 80;
N2 = 80;

epsilon1R=10^-5;
epsilon2R=10^-4;

Beta_1 = 2^(N1/M1)-1;
Beta_2 = 2^(N2/M1)-1;

delta_1 = sqrt(1/(2*pi*(2^(2*N1/M1)-1)));
delta_2 = sqrt(1/(2*pi*(2^(2*N2/M1)-1)));

nu_1 = Beta_1 - (1/(2*sqrt(M1)*delta_1));
mu_1 = Beta_1 + (1/(2*sqrt(M1)*delta_1));
nu_2 = Beta_2 - (1/(2*sqrt(M2)*delta_2));
mu_2 = Beta_2 + (1/(2*sqrt(M2)*delta_2));

% uniform random numbers
x = rand(N,1);

for u =1:length(RHO)
    
    % Received SINR
    % inverse of CDF
    Fgamma2_inv = @(t) (lambda.*rho(u).*alpha2(1).*log(1-t))./((lambda.*rho(u).*alpha1(1).*log(1-t))-2); 
    Fgamma12_inv = @(t) (lambda.*rho(u).*alpha2(1).*log(1-sqrt(t)))./((lambda.*rho(u).*alpha1(1).*log(1-sqrt(t)))-1);
    Fgamma11_inv = @(t) -lambda.*rho(u).*alpha1(1).*log(1-sqrt(t));
    
    % random numbers according to your CDF
    gamma2(u,:) = Fgamma2_inv(x);
    gamma12(u,:) = Fgamma12_inv(x);
    gamma11(u,:) = Fgamma11_inv(x);
   
    % Q-function
    x_2a(u,:) = log2 (1+gamma2(u,:)) - N2/M2;
    x_2b(u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma2(u,:)).^2)./M2);
    ratio2(u,:) = x_2a(u,:)./x_2b(u,:);
    q2(1,u,:) = qfunc(ratio2(u,:));
    
    x_12a(u,:) = log2 (1+gamma12(u,:)) - N2/M1;
    x_12b(u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma12(u,:)).^2)./M1);
    ratio12(u,:) = x_12a(u,:)./x_12b(u,:);
    q12(1,u,:) = qfunc(ratio12(u,:)); 
    
    x_11a(u,:) = log2 (1+gamma11(u,:)) - N1/M1;
    x_11b(u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma11(u,:)).^2)./M1);
    ratio11(u,:) = x_11a(u,:)./x_11b(u,:);
    q11(1,u,:) = qfunc(ratio11(u,:));
    
    % Received SINR
    % inverse of CDF
    Fgamma2_inv = @(t) (lambda.*rho(u).*alpha2(2).*log(1-t))./((lambda.*rho(u).*alpha1(2).*log(1-t))-2); 
    Fgamma12_inv = @(t) (lambda.*rho(u).*alpha2(2).*log(1-sqrt(t)))./((lambda.*rho(u).*alpha1(2).*log(1-sqrt(t)))-1);
    Fgamma11_inv = @(t) -lambda.*rho(u).*alpha1(2).*log(1-sqrt(t));
    
    % random numbers according to your CDF
    gamma2(u,:) = Fgamma2_inv(x);
    gamma12(u,:) = Fgamma12_inv(x);
    gamma11(u,:) = Fgamma11_inv(x);
   
    % Q-function
    x_2a(u,:) = log2 (1+gamma2(u,:)) - N2/M2;
    x_2b(u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma2(u,:)).^2)./M2);
    ratio2(u,:) = x_2a(u,:)./x_2b(u,:);
    q2(2,u,:) = qfunc(ratio2(u,:));
    
    x_12a(u,:) = log2 (1+gamma12(u,:)) - N2/M1;
    x_12b(u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma12(u,:)).^2)./M1);
    ratio12(u,:) = x_12a(u,:)./x_12b(u,:);
    q12(2,u,:) = qfunc(ratio12(u,:)); 
    
    x_11a(u,:) = log2 (1+gamma11(u,:)) - N1/M1;
    x_11b(u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma11(u,:)).^2)./M1);
    ratio11(u,:) = x_11a(u,:)./x_11b(u,:);
    q11(2,u,:) = qfunc(ratio11(u,:));
                

    % BLER from approximation Q + Riemann integral
    epsilon2_high_SNR(u) = (2 * Beta_2)/(lambda*rho(u)*(alpha2(1)-alpha1(1) * Beta_2));
    epsilon1_high_SNR(u) = (Beta_1/(lambda*rho(u)*alpha1(1)))^2 + (epsilon2R / 2)^2;
    
    epsilon2(u) = 1 - exp((-2 * Beta_2)/(lambda*rho(u)*(alpha2(1)-alpha1(1) * Beta_2)));
    epsilon1(u) = (1-exp(-Beta_1/(lambda*rho(u)*alpha1(1))))^2+(1-exp(-Beta_2/(lambda*rho(u)*(alpha2(1)-alpha1(1)*Beta_2))))^2;
    
end


figure (1)
semilogy(RHO, mean(q2(1,:,:),3),'b', 'linewidth', 1.5);
hold on; grid on;
semilogy(RHO, mean(q12(1,:,:),3)+(1-mean(q12(1,:,:),3)).*mean(q11(1,:,:),3),'r', 'linewidth', 1.5);



semilogy(RHO, epsilon2_high_SNR,'*b');
semilogy(RHO, epsilon1_high_SNR,'*r');

semilogy(RHO, epsilon2,'+b');
semilogy(RHO, epsilon1,'+r');

semilogy(RHO, mean(q2(2,:,:),3),'--b', 'linewidth', 1.5);
semilogy(RHO, mean(q12(2,:,:),3)+(1-mean(q12(2,:,:),3)).*mean(q11(2,:,:),3),'--r', 'linewidth', 1.5);
% axis([0 30 1e-5 1]);
title('BLER vs Transmit SNR');
xlabel('Transmit SNR (Rho in dB)');
ylabel('BLER');
dim = [0.2 0.5 0.01 0.01];
str = {'Solida line : \alpha_1=0.1','Dashed line : \alpha_1=0.2'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
% 
% legend('epsilon far1','epsilon near1','epsilon far2','epsilon near2',...
%        'epsilon far1 high SINR','epsilon near1 high SINR',...
%        'epsilon far2 appro','epsilon near2 appro');


legend('epsilon far1','epsilon near1',...
       'epsilon far1 high SINR','epsilon near1 high SINR',...
       'epsilon far2 appro','epsilon near2 appro');

figure (2)


semilogy(RHO, mean(q2(1,:,:),3),'b', 'linewidth', 1.5);
hold on; grid on;
semilogy(RHO, mean(q12(1,:,:),3)+(1-mean(q12(1,:,:),3)).*mean(q11(1,:,:),3),'r', 'linewidth', 1.5);
semilogy(RHO, mean(q2(2,:,:),3),'--b', 'linewidth', 1.5);
semilogy(RHO, mean(q12(2,:,:),3)+(1-mean(q12(2,:,:),3)).*mean(q11(2,:,:),3),'--r', 'linewidth', 1.5);

axis([0 0 1e-30 1]);
title('BLER vs Transmit SNR');
xlabel('Transmit SNR (Rho in dB)');
ylabel('BLER');
dim = [0.2 0.5 0.01 0.01];
str = {'Solida line : \alpha_1=0.1','Dashed line : \alpha_1=0.2'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
% 
% legend('epsilon far1','epsilon near1','epsilon far2','epsilon near2',...
%        'epsilon far1 high SINR','epsilon near1 high SINR',...
%        'epsilon far2 appro','epsilon near2 appro');


legend('epsilon far1','epsilon near1');

   

toc









