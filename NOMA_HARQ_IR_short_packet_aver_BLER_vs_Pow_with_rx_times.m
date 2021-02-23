clc; clear variables; close all;
tic
% Basic setting
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
M1 = [200 400 600];
M2 = [200 400 600];
N1 = 256;
N2 = 256;


% Transmit power in dBm
Pt = 0:2:30;                
% Pt = -113:2:-74;
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);

% Power allocaitn coefficient
alpha1 = 0.05;
alpha2 = 1 - alpha1;

gamma2 = zeros(length(Pt),N);
gamma12 = zeros(length(Pt),N);
gamma11 = zeros(length(Pt),N);

for mm = 1:length(M1)
    parfor u=1:length(Pt)
        rho = (10^-3)*db2pow(Pt(u))/ no;
        % Q-function w/ NO approximation
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
        
        % Common term
        w1 = 2^(N1/M1(mm))-1;
        w2 = 2^(N2/M2(mm))-1;

        Xi_1 = sqrt(1/(2*pi*(2^(2*N1/M1(mm))-1)));
        Xi_2 = sqrt(1/(2*pi*(2^(2*N2/M2(mm))-1)));

        nu_1 = w1 - (1/(2*sqrt(M1(mm))*Xi_1));
        tau_1 = w1 + (1/(2*sqrt(M1(mm))*Xi_1));
        nu_2 = w2 - (1/(2*sqrt(M2(mm))*Xi_2));
        tau_2 = w2 + (1/(2*sqrt(M2(mm))*Xi_2));
        
        % BLER from linerization Q
        first_term_22 = ei(-alpha2/(lambda2*rho*alpha1*(alpha2-alpha1*nu_2)))-...
        ei(-alpha2/(lambda2*rho*alpha1*(alpha2-alpha1*tau_2)));

        first_term_12 = ei(-alpha2/(lambda1*rho*alpha1*(alpha2-alpha1*nu_2)))-...
        ei(-alpha2/(lambda1*rho*alpha1*(alpha2-alpha1*tau_2)));

        second_term_22 = exp(-nu_2/(lambda2*rho*(alpha2-alpha1*nu_2)))*(alpha2-alpha1*nu_2) -...
        exp(-tau_2/(lambda2*rho*(alpha2-alpha1*tau_2)))*(alpha2-alpha1*tau_2);

        second_term_12 = exp(-nu_2/(lambda1*rho*(alpha2-alpha1*nu_2)))*(alpha2-alpha1*nu_2) - ...
        exp(-tau_2/(lambda1*rho*(alpha2-alpha1*tau_2)))*(alpha2-alpha1*tau_2); 

        Epsilon2(mm,u) = 1 - ((alpha2*Xi_2*sqrt(M2(mm))*exp(1/(lambda2*rho*alpha1)))/(lambda2*rho*(alpha1^2)))* first_term_22...
        - ((Xi_2*sqrt(M2(mm)))/alpha1)*second_term_22;

        Epsilon12(mm,u) = 1 - ((alpha2*Xi_2*sqrt(M2(mm))*exp(1/(lambda1*rho*alpha1)))/(lambda1*rho*(alpha1^2)))* first_term_12...
        - ((Xi_2*sqrt(M2(mm)))/alpha1)*second_term_12;


        Epsilon11(mm,u) = 1 + Xi_1*sqrt(M1(mm))*lambda1*alpha1*rho*(exp(-tau_1/(lambda1*alpha1*rho))-exp(-nu_1/(lambda1*alpha1*rho)));
        Epsilon1(mm,u) = Epsilon12(mm,u) + Epsilon11(mm,u);   
        
        %% BLER from approximation Q + Riemann integral + High SNR approximation
        epsilon2_high_SNR(mm,u) = w2/(lambda2*rho*(alpha2-alpha1*w2));
        epsilon12_high_SNR(mm,u) = w2/(lambda1*rho*(alpha2-alpha1*w1));
        epsilon11_high_SNR(mm,u) = w1/(lambda1*rho*alpha1);
        epsilon1_high_SNR(mm,u) = epsilon12_high_SNR(mm,u) + epsilon11_high_SNR(mm,u);


    end
end




% Plot

figure (1)
semilogy(Pt, mean(q2(1,:,:),3),'r');
hold on; grid on;
semilogy(Pt, mean(q12(1,:,:),3)+(1-mean(q12(1,:,:),3)).*mean(q11(1,:,:),3),'b');

semilogy(Pt, mean(q2(2,:,:),3),'--r');
semilogy(Pt, mean(q12(2,:,:),3)+(1-mean(q12(2,:,:),3)).*mean(q11(2,:,:),3),'--b');

semilogy(Pt, mean(q2(3,:,:),3),'-.r');
semilogy(Pt, mean(q12(3,:,:),3)+(1-mean(q12(3,:,:),3)).*mean(q11(3,:,:),3),'-.b'    );

semilogy(Pt, Epsilon2(1,:),'*r');
semilogy(Pt, Epsilon1(1,:),'*b');
semilogy(Pt, epsilon2_high_SNR(1,:),'or');
semilogy(Pt, epsilon1_high_SNR(1,:),'ob');

semilogy(Pt, Epsilon2(2,:),'*r');
semilogy(Pt, Epsilon1(2,:),'*b');
semilogy(Pt, epsilon2_high_SNR(2,:),'or');
semilogy(Pt, epsilon1_high_SNR(2,:),'ob');

semilogy(Pt, Epsilon2(3,:),'*r');
semilogy(Pt, Epsilon1(3,:),'*b');
semilogy(Pt, epsilon2_high_SNR(3,:),'or');
semilogy(Pt, epsilon1_high_SNR(3,:),'ob');

% title('BLER vs Transmit Power');
xlabel('Transmit Power (dBm)');
ylabel('Average BLER');

x = [0.3 0.5];
y = [0.6 0.5];
an = annotation('textarrow',x,y,'String','T=1');
an.FontName = 'Times New Roman';

x = [0.7 0.5];
y = [0.8 0.5];
an2 = annotation('textarrow',x,y,'String','T=2');
an2.FontName = 'Times New Roman';

x = [0.1 0.5];
y = [0.2 0.5];
an3 = annotation('textarrow',x,y,'String','T=3');
an3.FontName = 'Times New Roman';

cn1 = annotation('ellipse',[.84 .68 .05 .05])
cn2 = annotation('ellipse',[.84 .68 .05 .05])
cn3 = annotation('ellipse',[.84 .68 .05 .05])

legend('Far user','Near user');

set(gca, 'FontName', 'Times New Roman');
toc