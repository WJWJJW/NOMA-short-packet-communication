clc; clear variables; close all;
tic
%% Basic setting
N = 1e6;
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

% AWGN noise
% BW = 10^6;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

% Parameter setting
M1 = [100 200 300];
M2 = [100 200 300];
N1 = 80;
N2 = 80;

w1 = 2.^(N1./M1)-1;
w2 = 2.^(N2./M2)-1;

% Transmit power in dBm
Pt = 0:2:30;                
% Pt = -113:2:-74;
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);

% Power allocaitn coefficient
alpha1 = 0.2;
alpha2 = 1 - alpha1;
%%
for mm = 1:length(M1)
    for u=1:length(Pt)
        rho = (10^-3)*db2pow(Pt(u))/ no;
        %% Q-function w/ NO approximation
        % Received SINR
        gamma2(u,:) = (alpha2*pt(u).*abs(h2).^2)./(alpha1*pt(u).*abs(h2).^2+no);
        gamma12(u,:) = (alpha2*pt(u).*abs(h1).^2)./(alpha1*pt(u).*abs(h1).^2+no);
        gamma11(u,:) = (alpha1*pt(u).*abs(h1).^2)./no;

        % Q-function argument
        x_2a(mm,u,:) = log2 (1+gamma2(u,:)).*M2(mm) - N2;
        x_2b(mm,u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma2(u,:)).^2) .*M2(mm));
        ratio2(mm,u,:) = x_2a(mm,u,:)./x_2b(mm,u,:);

        x_12a(mm,u,:) = log2 (1+gamma12(u,:)).*M1(mm) - N2;
        x_12b(mm,u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma12(u,:)).^2) .*M1(mm));
        ratio12(mm,u,:) = x_12a(mm,u,:)./x_12b(mm,u,:);

        x_11a(mm,u,:) = log2 (1+gamma11(u,:)).*M1(mm) - N1;
        x_11b(mm,u,:) = log2 (exp(1))*sqrt((1-1./(1+gamma11(u,:)).^2) .*M1(mm));
        ratio11(mm,u,:) = x_11a(mm,u,:)./x_11b(mm,u,:);

        % Q-function
        q2(mm,u,:) = qfunc(ratio2(mm,u,:));
        q12(mm,u,:) = qfunc(ratio12(mm,u,:)); 
        q11(mm,u,:) = qfunc(ratio11(mm,u,:));
        
%         %% BLER from approximation Q + Riemann integral + High SNR approximation
%         epsilon2_high_SNR(mm,u) = w2(mm)/(lambda2*rho*(alpha2-alpha1*w2(mm)));
%         epsilon12_high_SNR(mm,u) = w2(mm)/(lambda1*rho*(alpha2-alpha1*w1(mm)));
%         epsilon11_high_SNR(mm,u) = w1(mm)/(lambda1*rho*alpha1);
%         epsilon1_high_SNR(mm,u) = epsilon12_high_SNR(mm,u) + epsilon11_high_SNR(mm,u);


    end
end




%% Plot

figure (1)
semilogy(Pt, mean(q2(1,:,:),3),'r', 'linewidth', 1.5);
hold on; grid on;
semilogy(Pt, mean(q12(1,:,:),3)+(1-mean(q12(1,:,:),3)).*mean(q11(1,:,:),3),'b', 'linewidth', 1.5);
% semilogy(Pt, epsilon2_high_SNR(1,:),'or');
% semilogy(Pt, epsilon1_high_SNR(1,:),'ob');



semilogy(Pt, mean(q2(2,:,:),3),'--r', 'linewidth', 1.5);
semilogy(Pt, mean(q12(2,:,:),3)+(1-mean(q12(2,:,:),3)).*mean(q11(2,:,:),3),'--b', 'linewidth', 1.5);
% semilogy(Pt, epsilon2_high_SNR(2,:),'or');
% semilogy(Pt, epsilon1_high_SNR(2,:),'ob');

semilogy(Pt, mean(q2(3,:,:),3),'-.r', 'linewidth', 1.5);
semilogy(Pt, mean(q12(3,:,:),3)+(1-mean(q12(3,:,:),3)).*mean(q11(3,:,:),3),'-.b', 'linewidth', 1.5);
% semilogy(Pt, epsilon2_high_SNR(3,:),'or');
% semilogy(Pt, epsilon1_high_SNR(3,:),'ob');


title('BLER vs Transmit Power');
xlabel('Transmit Power (in dBm)');
ylabel('BLER');



legend('Far user','Near user');


toc