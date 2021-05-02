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
tx_times = 1:3;
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
epsilon1 = zeros(length(tx_times), length(tx_times)+1, length(Pt));

cc_epsilon2 = zeros(length(tx_times), length(tx_times)+1, length(Pt));
cc_epsilon1 = zeros(length(tx_times), length(tx_times)+1, length(Pt));

epsilon2(:,1,:) = 1;
epsilon1(:,1,:) = 1;

cc_epsilon2(:,1,:) = 1;
cc_epsilon1(:,1,:) = 1;

T_cu_2 = zeros(length(tx_times), length(Pt),length(Delay_related));
T_cu_1 = zeros(length(tx_times), length(Pt),length(Delay_related));

N_expected_2 = zeros(length(tx_times), length(Pt));
N_expected_1 = zeros(length(tx_times), length(Pt));

throughput_2 = zeros(length(tx_times), length(Pt),length(Delay_related));
throughput_1 = zeros(length(tx_times), length(Pt),length(Delay_related));

cc_T_cu_2 = zeros(length(tx_times), length(Pt),length(Delay_related));
cc_T_cu_1 = zeros(length(tx_times), length(Pt),length(Delay_related));

cc_N_expected_2 = zeros(length(tx_times), length(Pt));
cc_N_expected_1 = zeros(length(tx_times), length(Pt));

cc_throughput_2 = zeros(length(tx_times), length(Pt),length(Delay_related));
cc_throughput_1 = zeros(length(tx_times), length(Pt),length(Delay_related));

for u=1:length(Pt)
    rho = pt(u)/ no;
    for tt = 1:length(tx_times)
        cc_gamma2 = 0;
        cc_gamma12 = 0;
        cc_gamma11 = 0;
        sub_M = M / tt;
        inc_M = M / tt;
        for mm = 1:tt
            % Q-function w/ NO approximation
            % Received SINR (HARQ-IR)
            gamma2 = (alpha2*pt(u).*abs(h2).^2)./(alpha1*pt(u).*abs(h2).^2+no);
            gamma12 = (alpha2*pt(u).*abs(h1).^2)./(alpha1*pt(u).*abs(h1).^2+no);
            gamma11 = (alpha1*pt(u).*abs(h1).^2)./no;
            
            % Received SINR (HARQ-CC)
            cc_gamma2 = cc_gamma2 + gamma2;
            cc_gamma12 = cc_gamma12 + gamma12;
            cc_gamma11 = cc_gamma11 + gamma11;

            % Q-function argument (HARQ-IR)
            x_2a = log2 (1+gamma2).*sub_M - N2;
            x_2b = log2 (exp(1))*sqrt((1-1./(1+gamma2).^2) .*sub_M);
            ratio2 = x_2a./x_2b;

            x_12a = log2 (1+gamma12).*sub_M - N2;
            x_12b = log2 (exp(1))*sqrt((1-1./(1+gamma12).^2) .*sub_M);
            ratio12 = x_12a./x_12b;

            x_11a = log2 (1+gamma11).*sub_M - N1;
            x_11b = log2 (exp(1))*sqrt((1-1./(1+gamma11).^2) .*sub_M);
            ratio11 = x_11a./x_11b;

            % Q-function (HARQ-IR)
            q2 = qfunc(ratio2);
            q12 = qfunc(ratio12); 
            q11 = qfunc(ratio11);
            
            % Q-function argument (HARQ-CC)
            x_2a = log2 (1+cc_gamma2).*M - N2;
            x_2b = log2 (exp(1))*sqrt((1-1./(1+cc_gamma2).^2) .*M);
            ratio2 = x_2a./x_2b;

            x_12a = log2 (1+cc_gamma12).*M - N2;
            x_12b = log2 (exp(1))*sqrt((1-1./(1+cc_gamma12).^2) .*M);
            ratio12 = x_12a./x_12b;

            x_11a = log2 (1+cc_gamma11).*M - N1;
            x_11b = log2 (exp(1))*sqrt((1-1./(1+cc_gamma11).^2) .*M);
            ratio11 = x_11a./x_11b;
            
            % Q-function (HARQ-CC)
            cc_q2 = qfunc(ratio2);
            cc_q12 = qfunc(ratio12); 
            cc_q11 = qfunc(ratio11);
            
            % Average BLER (HARQ-IR)
            epsilon2(tt,mm+1,u) = mean(q2);
            epsilon1(tt,mm+1,u) = mean(q12)+(1-mean(q12))*mean(q11);
            
            % Average BLER (HARQ-CC)
            cc_epsilon2(tt,mm+1,u) = mean(cc_q2);
            cc_epsilon1(tt,mm+1,u) = mean(cc_q12)+(1-mean(cc_q12))*mean(cc_q11);
                
            sub_M = sub_M + inc_M;
        end
        if tt == 1
            T_cu_2(tt,u,:) = M*epsilon2(tt,1,u);
            T_cu_1(tt,u,:) = M*epsilon1(tt,1,u);
            cc_T_cu_2(tt,u,:) = M*cc_epsilon2(tt,1,u);
            cc_T_cu_1(tt,u,:) = M*cc_epsilon1(tt,1,u);
        else
            T_cu_2(tt,u,:) = sum(inc_M*epsilon2(tt,1:tt,u))+...
                                    Delay.*sum(epsilon2(tt,1:tt,u));
            T_cu_1(tt,u,:) = sum(inc_M*epsilon1(tt,1:tt,u))+...
                                    Delay.*sum(epsilon1(tt,1:tt,u));  
            cc_T_cu_2(tt,u,:) = sum(M*cc_epsilon2(tt,1:tt,u))+...
                                    Delay.*sum(cc_epsilon2(tt,1:tt,u));
            cc_T_cu_1(tt,u,:) = sum(M*cc_epsilon1(tt,1:tt,u))+...
                                    Delay.*sum(cc_epsilon1(tt,1:tt,u));  
        end
        N_expected_2(tt,u) = N2*(1-epsilon2(tt,mm+1,u));
        N_expected_1(tt,u) = N1*(1-epsilon1(tt,mm+1,u));     
        cc_N_expected_2(tt,u) = N2*(1-cc_epsilon2(tt,mm+1,u));
        cc_N_expected_1(tt,u) = N1*(1-cc_epsilon1(tt,mm+1,u));  
    end
    % Calculate throughput
    throughput_2(:,u,:) = N_expected_2(:,u) ./ T_cu_2(:,u,:);
    throughput_1(:,u,:) = N_expected_1(:,u) ./ T_cu_1(:,u,:);
    cc_throughput_2(:,u,:) = cc_N_expected_2(:,u) ./ cc_T_cu_2(:,u,:);
    cc_throughput_1(:,u,:) = cc_N_expected_1(:,u) ./ cc_T_cu_1(:,u,:);
end

figure (1)
plot(Pt, throughput_2(1,:,1),'b');
hold on; grid on;
plot(Pt, throughput_2(3,:,1),'r');
plot(Pt, cc_throughput_2(3,:,1),'g');

plot(Pt, throughput_2(1,:,2),'ob');
plot(Pt, throughput_2(1,:,3),'*b');


plot(Pt, throughput_2(3,:,2),'or');
plot(Pt, throughput_2(3,:,3),'*r');


plot(Pt, cc_throughput_2(3,:,2),'og');
plot(Pt, cc_throughput_2(3,:,3),'*g');

xlabel('Transmitted power (dBm)');
ylabel('Throughput (bpcu)');
legend('One-shot','HARQ-IR','HARQ-CC');
set(gca, 'FontName', 'Times New Roman');



figure (2)
plot(Pt, throughput_1(1,:,1),'b');
hold on; grid on;
plot(Pt, throughput_1(3,:,1),'r');
plot(Pt, cc_throughput_1(3,:,1),'g');


plot(Pt, throughput_1(1,:,2),'ob');
plot(Pt, throughput_1(1,:,3),'*b');


plot(Pt, throughput_1(3,:,2),'or');
plot(Pt, throughput_1(3,:,3),'*r');


plot(Pt, cc_throughput_1(3,:,2),'og');
plot(Pt, cc_throughput_1(3,:,3),'*g');

xlabel('Transmitted power (dBm)');
ylabel('Throughput (bpcu)');
legend('One-shot','HARQ-IR','HARQ-CC');
set(gca, 'FontName', 'Times New Roman');


figure (3)
plot(Pt, T_cu_2(1,:,1),'b');
hold on; grid on;
plot(Pt, T_cu_2(3,:,1),'r');
plot(Pt, cc_T_cu_2(3,:,1),'g');


plot(Pt, T_cu_2(1,:,2),'ob');
plot(Pt, T_cu_2(1,:,3),'*b');


plot(Pt, T_cu_2(3,:,2),'or');
plot(Pt, T_cu_2(3,:,3),'*r');


plot(Pt, cc_T_cu_2(3,:,2),'og');
plot(Pt, cc_T_cu_2(3,:,3),'*g');

xlabel('Transmitted power (dBm)');
ylabel('Expected delay (cu)');
legend('One-shot','HARQ-IR','HARQ-CC');
set(gca, 'FontName', 'Times New Roman');


figure (4)
plot(Pt, T_cu_1(1,:,1),'b');
hold on; grid on;
plot(Pt, T_cu_1(3,:,1),'r');
plot(Pt, cc_T_cu_1(3,:,1),'g');


plot(Pt, T_cu_1(1,:,2),'ob');
plot(Pt, T_cu_1(1,:,3),'*b');


plot(Pt, T_cu_1(3,:,2),'or');
plot(Pt, T_cu_1(3,:,3),'*r');


plot(Pt, cc_T_cu_1(3,:,2),'og');
plot(Pt, cc_T_cu_1(3,:,3),'*g');

xlabel('Transmitted power (dBm)');
ylabel('Expected delay (cu)');
legend('One-shot','HARQ-IR','HARQ-CC');
set(gca, 'FontName', 'Times New Roman');

figure (5)
plot(Pt,reshape(epsilon2(1,2,:),1,length(Pt)),'b');
hold on; grid on;
plot(Pt,reshape(epsilon2(3,2,:),1,length(Pt)),'r');
plot(Pt,reshape(cc_epsilon2(3,2,:),1,length(Pt)),'--g');


plot(Pt,reshape(epsilon2(3,3,:),1,length(Pt)),'-or');
plot(Pt,reshape(epsilon2(3,4,:),1,length(Pt)),'*r');


plot(Pt,reshape(cc_epsilon2(3,3,:),1,length(Pt)),'-og');
plot(Pt,reshape(cc_epsilon2(3,4,:),1,length(Pt)),'-*g');


xlabel('Transmitted power (dBm)');
ylabel('BLER');
legend('One-shot','HARQ-IR','HARQ-CC');
set(gca, 'FontName', 'Times New Roman');

figure (6)
plot(Pt,reshape(epsilon1(1,2,:),1,length(Pt)),'b');
hold on; grid on;
plot(Pt,reshape(epsilon1(3,2,:),1,length(Pt)),'r');
plot(Pt,reshape(cc_epsilon1(3,2,:),1,length(Pt)),'--g');

plot(Pt,reshape(epsilon1(3,3,:),1,length(Pt)),'-or');
plot(Pt,reshape(epsilon1(3,4,:),1,length(Pt)),'*r');

plot(Pt,reshape(cc_epsilon1(3,3,:),1,length(Pt)),'-og');
plot(Pt,reshape(cc_epsilon1(3,4,:),1,length(Pt)),'-*g');


xlabel('Transmitted power (dBm)');
ylabel('BLER');
legend('One-shot','HARQ-IR','HARQ-CC');
set(gca, 'FontName', 'Times New Roman');