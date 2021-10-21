% This script analysis relation between power and BLER
% fixed blocklength and total transmitted power for OMA user

% This script analysis relation between power and BLER
% fixed blocklength and total transmitted power

clc; clear variables; close all;
tic
%% Basic setting
% # of channel tap
N = 1e7;
% Path loss exponent
eta = 4;
% User distance
d1 = 100;
d2 = 200;
% Rayleigh fading channel
h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));
% Channel mean
lambda1 = mean(abs(h1).^2);
lambda2 = mean(abs(h2).^2);

% Transmit power in dBm
Pt = 0:2:30;                
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);

% AWGN noise
% BW = 10^6;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

% # of channel use (blocklength)
M1 = 300;
M2 = 300;
% M = 150;

% # of information bits
N1 = 256;
N2 = 256;

%DoF coefficient
beta1 = 1/2;
beta2 = 1 - beta1;

% Common term
w1o = 2^(N1/(M1*beta1))-1;
w2o = 2^(N2/(M2*beta2))-1;

Xi_1 = sqrt(1/(2*pi*(2^(2*N1/(M1*beta1))-1)));
Xi_2 = sqrt(1/(2*pi*(2^(2*N2/(M2*beta2))-1)));

nu_1 = w1o - (1/(2*sqrt(M1)*Xi_1));
tau_1 = w1o + (1/(2*sqrt(M1)*Xi_1));
nu_2 = w2o - (1/(2*sqrt(M2)*Xi_2));
tau_2 = w2o + (1/(2*sqrt(M2)*Xi_2));

%%

for u=1:length(Pt)
    rho = pt(u)/ no;
    %% Q-function w/ NO approximation
    % Received SINR
    gamma2(u,:) = (pt(u).*abs(h2).^2)./(beta2*no);
    gamma1(u,:) = (pt(u).*abs(h1).^2)./(beta1*no);
    
    % Q-function argument
    x_2a(u,:) = sqrt(M2).*(log2 (1+gamma2(u,:)) - N2/(M2*beta2));
    x_2b(u,:) = log2 (exp(1))*sqrt(1-1./(1+gamma2(u,:)).^2);
    ratio2(u,:) = x_2a(u,:)./x_2b(u,:);

    x_1a(u,:) = sqrt(M1).*(log2 (1+gamma1(u,:)) - N1/(M1*beta1));
    x_1b(u,:) = log2 (exp(1))*sqrt(1-1./(1+gamma1(u,:)).^2);
    ratio1(u,:) = x_1a(u,:)./x_1b(u,:);
    
    % Q-function
    q2(u,:) = qfunc(ratio2(u,:)); 
    q1(u,:) = qfunc(ratio1(u,:)); 
    

    %% BLER from approximation Q + Riemann integral
    epsilon2(u) = 1 - exp(-(w2o*beta2)/(lambda2*rho));
    epsilon1(u) = 1 - exp(-(w1o*beta1)/(lambda1*rho));
    
    
    %% BLER from approximation Q + Riemann integral + High SNR approximation
    epsilon2_high_SNR(u) = (w2o*beta2)/(lambda2*rho);
    epsilon1_high_SNR(u) = (w1o*beta1)/(lambda1*rho);;
    
end

%% Plot

figure (1)
semilogy(Pt, mean(q2,2),'r', 'linewidth', 1.5);
hold on; grid on;
semilogy(Pt, mean(q1,2),'b', 'linewidth', 1.5);

semilogy(Pt, epsilon2,'or');
semilogy(Pt, epsilon1, 'ob');


semilogy(Pt, epsilon2_high_SNR,'*r');
semilogy(Pt, epsilon1_high_SNR, '*b');


% axis([0 40 1e-6 1]);
title('BLER vs Transmitted power');
xlabel('Transmitted power (in dBm)');
ylabel('BLER');
dim = [0.15 0.5 0.005 0.01];
str = {'Line : No approximation',...
       'Circle :Q-function linearization + ','Riemann integral',...
       'Star : Q-function linearization + ','Riemann integral + ','High SNR approximation'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
legend('BLER2','BLER1');

toc

