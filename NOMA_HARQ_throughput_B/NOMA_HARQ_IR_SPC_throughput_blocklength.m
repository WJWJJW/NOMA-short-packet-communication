clc; clear variables; close all;

% Basic setting
% # of channel tap
N = 1e6;

% Path loss exponent
eta = 4;
% User distance
d1 = 100;
d2 = 200;

% Transmit power in dBm
Pt = 30;                
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);

% AWGN noise
% BW = 10^6;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -80;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt / no;

% User distance
h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));


alpha1 = 0.2;
alpha2 = 1-alpha1;

% # of channel use (blocklength)
tx_times = 1:2;

% # of information bits
N1 = 256;
N2 = 256;

Delay_related = [0 0.05 0.1];

M = 200:20:1000;

epsilon2 = zeros(length(tx_times), length(tx_times)+1, length(M));
epsilon1 = zeros(length(tx_times), length(tx_times)+1, length(M));

epsilon2(:,1,:) = 1;
epsilon1(:,1,:) = 1;

T_cu_2 = zeros(length(tx_times), length(M),length(Delay_related));
T_cu_1 = zeros(length(tx_times), length(M),length(Delay_related));

N_expected_2 = zeros(length(tx_times), length(M));
N_expected_1 = zeros(length(tx_times), length(M));

throughput_2 = zeros(length(tx_times), length(M),length(Delay_related));
throughput_1 = zeros(length(tx_times), length(M),length(Delay_related));

for u=1:length(M)
    % Received SINR
    gamma2 = (alpha2*pt.*abs(h2).^2)./(alpha1*pt.*abs(h2).^2+no);
    gamma12 = (alpha2*pt.*abs(h1).^2)./(alpha1*pt.*abs(h1).^2+no);
    gamma11 = (alpha1*pt.*abs(h1).^2)./no;
    
    
    MM = M(u);
    % Delay
    Delay = Delay_related.*400;
        for tt = 1:length(tx_times)
            sub_M = MM;
            inc_M = MM;
            for mm = 1:tt
                % Q-function w/ NO approximation
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
                T_cu_2(tt,u,:) = MM*epsilon2(tt,1,u);
                T_cu_1(tt,u,:) = MM*epsilon1(tt,1,u);
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
plot(M, throughput_2(1,:,1),'b');
hold on; grid on;
plot(M, throughput_2(1,:,2),'ob');
plot(M, throughput_2(1,:,3),'*b');

plot(M, throughput_2(2,:,1),'r');
plot(M, throughput_2(2,:,2),'-or');
plot(M, throughput_2(2,:,3),'-*r');

xlabel('Blocklength (Channel uses)');
ylabel('Throughput (bpcu)');

set(gca, 'FontName', 'Times New Roman');

figure (2)
plot(M, throughput_1(1,:,1),'b');
hold on; grid on;
plot(M, throughput_1(1,:,2),'ob');
plot(M, throughput_1(1,:,3),'*b');

plot(M, throughput_1(2,:,1),'r');
plot(M, throughput_1(2,:,2),'-or');
plot(M, throughput_1(2,:,3),'-*r');

xlabel('Blocklength (Channel uses)');
ylabel('Throughput (bpcu)');

set(gca, 'FontName', 'Times New Roman');




