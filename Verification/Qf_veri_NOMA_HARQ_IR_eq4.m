clc; clear variables; close all;
tic
% Basic setting
N = 1e6;
eta = 4;
% User distance
d1 = 400;
d2 = 900;

% Reliability constraint
eplsion1R = 10^-5;
eplsion2R = 10^-4;

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
N1 = 80;
N2 = 80;

% Transmit power in dBm
Pt = 40;                
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);

% Power allocaitn coefficient
alpha1 = 0.1:0.1:0.5;
alpha2 = 1 - alpha1;

kk = randi([1 N],1);

for aa=1:length(alpha1)
    syms M1
    syms M2
    
    % Q-function w/ NO approximation
    % Received SINR
    gamma2(aa) = (alpha2(aa)*pt*abs(h2(kk))^2) / (alpha1(aa)*pt*abs(h2(kk))^2+no);
    gamma12(aa) = (alpha2(aa)*pt*abs(h1(kk))^2) / (alpha1(aa)*pt*abs(h1(kk))^2+no);
    gamma11(aa) = (alpha1(aa)*pt*abs(h1(kk))^2) / no;
    
    % Q-function argument
    x_2a(aa) = log2 (1+gamma2(aa))*M2 - N2;
    x_2b(aa) = log2 (exp(1))*sqrt((1-1/(1+gamma2(aa))^2) * M2);
    ratio2(aa) = x_2a(aa) / x_2b(aa);
    
    x_12a(aa) = log2 (1+gamma12(aa))*M1 - N2;
    x_12b(aa) = log2 (exp(1))*sqrt((1-1/(1+gamma12(aa))^2) * M1);
    ratio12(aa) = x_12a(aa) / x_12b(aa);
    
    x_11a(aa) = log2 (1+gamma11(aa))*M1 - N1;
    x_11b(aa) = log2 (exp(1))*sqrt((1-1/(1+gamma11(aa))^2) * M1);
    ratio11(aa) = x_11a(aa) / x_11b(aa);
    
    % Q-function
    q2(aa) = 0.5 * erfc(ratio2(aa)/sqrt(2));
    q12(aa) = 0.5 * erfc(ratio12(aa)/sqrt(2));
    q11(aa) = 0.5 * erfc(ratio11(aa)/sqrt(2));
    
        
    % Required blocklength solving
    eqn = eplsion1R ==  q12(aa)+(1-q12(aa))*q11(aa);
    MM1(aa) = double(vpasolve(eqn,M1));
    
    eqn2 = eplsion2R ==  q2(aa);
    MM2(aa) = double(vpasolve(eqn2,M2));
    
    
end    


%% Plot
figure(1)
plot(alpha1, MM2, 'linewidth', 1.5);
hold on; grid on;
plot(alpha1, MM1,'-o', 'linewidth', 1.5);

% axis([0.1 0.5 -2000 10000])
xlabel('alpha1');
ylabel('M');
title('Required blocklength vs alpha1');
legend('M2','M1');
   

toc