clc; clear variables; close all;

% Basic setting
% # of channel tap
N = 1e6;
% Path loss exponent
eta = 4;
% User distance
d1 = 100;
d2 = 200;
% Rayleigh fading channel
h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));
% Channel mean
% lambda1 = mean(abs(h1).^2);
% lambda2 = mean(abs(h2).^2);

% Transmit power in dBm
Pt = 0:2:30;                
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);

% AWGN noise
% BW = 10^6;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -80;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

% # of channel use (blocklength)
tx_times = 1:2;
M = 400;

% # of information bits
N1 = 256;
N2 = 256;

Delay_related = [0 0.05 0.1];
Delay = Delay_related.*M;

% Power allocaitn coefficient
alpha1 = 0.2;
alpha2 = 1 - alpha1;

epsilon2 = zeros(length(tx_times), length(tx_times)+1, length(Pt));
epsilon12 = zeros(length(tx_times), length(tx_times)+1, length(Pt));
epsilon11 = zeros(length(tx_times), length(tx_times)+1, length(Pt));
epsilon1 = zeros(length(tx_times), length(tx_times)+1, length(Pt));

epsilon2(:,1,:) = 1;
epsilon12(:,1,:) = 1;
epsilon11(:,1,:) = 1;
epsilon1(:,1,:) = 1;

T_cu_2 = zeros(length(tx_times), length(Pt),length(Delay_related));
T_cu_1 = zeros(length(tx_times), length(Pt),length(Delay_related));

N_expected_2 = zeros(length(tx_times), length(Pt));
N_expected_1 = zeros(length(tx_times), length(Pt));

throughput_2 = zeros(length(tx_times), length(Pt),length(Delay_related));
throughput_1 = zeros(length(tx_times), length(Pt),length(Delay_related));

for u=1:length(Pt)
    rho = pt(u)/ no;
    for tt = 1:length(tx_times)
        sub_M = M / tt;
        inc_M = M / tt;
        for mm = 1:tt
            % Common term
%             w1 = 2^(N1 / sub_M)-1;
%             w2 = 2^(N2 / sub_M)-1;
            % Calculate BLER
%             epsilon2(tt,mm+1,u) = w2/(lambda2*rho*(alpha2-alpha1*w2));
%             epsilon12(tt,mm+1,u) = w2/(lambda1*rho*(alpha2-alpha1*w2));
%             epsilon11(tt,mm+1,u) = w1/(lambda1*rho*alpha1);
%             epsilon1(tt,mm+1,u) = ...
%                 epsilon12(tt,mm+1,u) + epsilon11(tt,mm+1,u);
            % Q-function w/ NO approximation
            % Received SINR
            gamma2 = (alpha2*pt(u).*abs(h2).^2)./(alpha1*pt(u).*abs(h2).^2+no);
            gamma12 = (alpha2*pt(u).*abs(h1).^2)./(alpha1*pt(u).*abs(h1).^2+no);
            gamma11 = (alpha1*pt(u).*abs(h1).^2)./no;

            % Q-function argument
            x_2a = log2 (1+gamma2).*sub_M - N2;
            x_2b = log2 (exp(1))*sqrt((1-1./(1+gamma2).^2) .*sub_M);
            ratio2 = x_2a./x_2b;

            x_12a = log2 (1+gamma12).*sub_M - N2;
            x_12b = log2 (exp(1))*sqrt((1-1./(1+gamma12).^2) .*sub_M);
            ratio12 = x_12a./x_12b;

            x_11a = log2 (1+gamma11).*sub_M - N1;
            x_11b = log2 (exp(1))*sqrt((1-1./(1+gamma11).^2) .*sub_M);
            ratio11 = x_11a./x_11b;

            % Q-function
            q2 = qfunc(ratio2);
            q12 = qfunc(ratio12); 
            q11 = qfunc(ratio11); 
            
            epsilon2(tt,mm+1,u) = mean(q2);
            epsilon1(tt,mm+1,u) = mean(q12)+(1-mean(q12))*mean(q11);
                
            sub_M = sub_M + inc_M;
        end
        if tt == 1
            T_cu_2(tt,u,:) = M*epsilon2(tt,1,u);
            T_cu_1(tt,u,:) = M*epsilon1(tt,1,u);
        else
            T_cu_2(tt,u,:) = sum(inc_M*epsilon2(tt,1:tt,u))+...
                                    Delay.*sum(epsilon2(tt,1:tt,u));
            T_cu_1(tt,u,:) = sum(inc_M*epsilon1(tt,1:tt,u))+...
                                    Delay.*sum(epsilon1(tt,1:tt,u));  
        end
        N_expected_2(tt,u) = N2*(1-epsilon2(tt,mm+1,u));
        N_expected_1(tt,u) = N1*(1-epsilon1(tt,mm+1,u));     
    end
    % Calculate throughput
    throughput_2(:,u,:) = N_expected_2(:,u) ./ T_cu_2(:,u,:);
    throughput_1(:,u,:) = N_expected_1(:,u) ./ T_cu_1(:,u,:);
end

figure (1)
plot(Pt, throughput_2(1,:,1),'b');
hold on; grid on;
plot(Pt, throughput_2(1,:,2),'ob');
plot(Pt, throughput_2(1,:,3),'*b');

plot(Pt, throughput_2(2,:,1),'r');
plot(Pt, throughput_2(2,:,2),'or');
plot(Pt, throughput_2(2,:,3),'*r');
 

xlabel('Transmitted power (dBm)');
ylabel('Throughput (bpcu)');
set(gca, 'FontName', 'Times New Roman');

figure (2)
plot(Pt, throughput_1(1,:,1),'b');
hold on; grid on;
plot(Pt, throughput_1(1,:,2),'ob');
plot(Pt, throughput_1(1,:,3),'*b');

plot(Pt, throughput_1(2,:,1),'r');
plot(Pt, throughput_1(2,:,2),'or');
plot(Pt, throughput_1(2,:,3),'*r');


xlabel('Transmitted power (dBm)');
ylabel('Throughput (bpcu)');

set(gca, 'FontName', 'Times New Roman');


figure (3)
plot(Pt, T_cu_2(1,:,1),'b');
hold on; grid on;
plot(Pt, T_cu_2(1,:,2),'ob');
plot(Pt, T_cu_2(1,:,3),'*b');

plot(Pt, T_cu_2(2,:,1),'r');
plot(Pt, T_cu_2(2,:,2),'or');
plot(Pt, T_cu_2(2,:,3),'*r');

xlabel('Transmitted power (dBm)');
ylabel('Expected delay (cu)');
set(gca, 'FontName', 'Times New Roman');


figure (4)
plot(Pt, T_cu_1(1,:,1),'b');
hold on; grid on;
plot(Pt, T_cu_1(1,:,2),'ob');
plot(Pt, T_cu_1(1,:,3),'*b');

plot(Pt, T_cu_1(2,:,1),'r');
plot(Pt, T_cu_1(2,:,2),'or');
plot(Pt, T_cu_1(2,:,3),'*r');

xlabel('Transmitted power (dBm)');
ylabel('Expected delay (cu)');
set(gca, 'FontName', 'Times New Roman');