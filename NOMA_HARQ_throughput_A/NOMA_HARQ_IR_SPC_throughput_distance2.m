clc; clear variables; close all;

% Basic setting
% # of channel tap
N = 1e6;
% Reliability constraint
eplsion1R = 10^-5;
eplsion2R = 10^-4;
% Path loss exponent
eta = 4;
% User distance
d1 = 100;
d2 = 150:5:500;

% Transmit power in dBm
Pt = 30;                
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);

% AWGN noise
% BW = 10^6;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt / no;

% find optimal delta
h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
lamda1 = mean(abs(h1).^2);
[delta] = delta_finder (lamda1*rho*eplsion1R);


% # of channel use (blocklength)
tx_times = 1:2;

% # of information bits
N1 = 256;
N2 = 256;

Delay_related = [0 0.3 0.4];


epsilon2 = zeros(length(tx_times), length(tx_times)+1, length(d2));
epsilon1 = zeros(length(tx_times), length(tx_times)+1, length(d2));

epsilon2(:,1,:) = 1;
epsilon1(:,1,:) = 1;

T_cu_2 = zeros(length(tx_times), length(d2),length(Delay_related));
T_cu_1 = zeros(length(tx_times), length(d2),length(Delay_related));

N_expected_2 = zeros(length(tx_times), length(d2));
N_expected_1 = zeros(length(tx_times), length(d2));

throughput_2 = zeros(length(tx_times), length(d2),length(Delay_related));
throughput_1 = zeros(length(tx_times), length(d2),length(Delay_related));

opt_a1 = zeros(1,length(d2));
opt_M = zeros(1,length(d2));

for dd=1:length(d2)
    % find optimal blocklength and power allocation coefficeint
    h2 = sqrt(1/2*d2(dd)^-eta)*(randn(1,N)+1i*randn(1,N));
    lamda2 = mean(abs(h2).^2);
    
    contraint = (lamda1*eplsion1R)/(lamda2*eplsion2R);
    theta = contraint * delta;
    
    % Adjust target BLER according to condition 2
    check = (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - (lamda2*eplsion2R.*rho + 4)*lamda2*eplsion2R > 0;
    
    switch (check)
        case 0
            eplsion2R_m = theta*eplsion2R;
        case 1
            eplsion2R_m = eplsion2R;
    end
    
    % optimal alpha1 and blocklength calculation
    opt_a1(dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
              4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m))) ...
              / (2*lamda2*eplsion2R_m*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m));
    opt_M(dd) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)))...
                         / (2*lamda2*eplsion2R_m)));
                     
    % Received SINR
    alpha1 = opt_a1(dd);
    alpha2 = 1-alpha1;
    gamma2 = (alpha2*pt.*abs(h2).^2)./(alpha1*pt.*abs(h2).^2+no);
    gamma12 = (alpha2*pt.*abs(h1).^2)./(alpha1*pt.*abs(h1).^2+no);
    gamma11 = (alpha1*pt.*abs(h1).^2)./no;
    
    
    M = opt_M(dd);
    % Delay
    Delay = Delay_related.*300;
        for tt = 1:length(tx_times)
            sub_M = M / tt;
            inc_M = M / tt;
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

                epsilon2(tt,mm+1,dd) = mean(q2);
                epsilon1(tt,mm+1,dd) = mean(q12)+(1-mean(q12))*mean(q11);

                sub_M = sub_M + inc_M;
            end
            if tt == 1
                T_cu_2(tt,dd,:) = M*epsilon2(tt,1,dd);
                T_cu_1(tt,dd,:) = M*epsilon1(tt,1,dd);
            else
                T_cu_2(tt,dd,:) = sum(inc_M*epsilon2(tt,1:tt,dd))+...
                                        Delay.*sum(epsilon2(tt,1:tt,dd));
                T_cu_1(tt,dd,:) = sum(inc_M*epsilon1(tt,1:tt,dd))+...
                                        Delay.*sum(epsilon1(tt,1:tt,dd));  
            end
            N_expected_2(tt,dd) = N2*(1-epsilon2(tt,mm+1,dd));
            N_expected_1(tt,dd) = N1*(1-epsilon1(tt,mm+1,dd));     
        end
        % Calculate throughput
        throughput_2(:,dd,:) = N_expected_2(:,dd) ./ T_cu_2(:,dd,:);
        throughput_1(:,dd,:) = N_expected_1(:,dd) ./ T_cu_1(:,dd,:);
end

syms D2
eqn_constraint2 = (D2^-eta*eplsion2R*rho + 2)*d1^-eta*eplsion1R - (D2^-eta.*eplsion2R.*rho + 4)*D2^-eta*eplsion2R == 0;
D2_sol2 = vpasolve(eqn_constraint2, D2);


figure (1)
plot(d2, throughput_2(1,:,1),'b');
hold on; grid on;
plot(d2, throughput_2(1,:,2),'ob');
plot(d2, throughput_2(1,:,3),'*b');

plot(d2, throughput_2(2,:,1),'r');
plot(d2, throughput_2(2,:,2),'-or');
plot(d2, throughput_2(2,:,3),'-*r');

xlabel('Distance (m)');
ylabel('Throughput (bpcu)');

xline(double(D2_sol2(2)),'-','Threshold1');
set(gca, 'FontName', 'Times New Roman');

figure (2)
plot(d2, throughput_1(1,:,1),'b');
hold on; grid on;
plot(d2, throughput_1(1,:,2),'ob');
plot(d2, throughput_1(1,:,3),'*b');

plot(d2, throughput_1(2,:,1),'r');
plot(d2, throughput_1(2,:,2),'-or');
plot(d2, throughput_1(2,:,3),'-*r');

xlabel('Distance (m)');
ylabel('Throughput (bpcu)');

xline(double(D2_sol2(2)),'-','Threshold1');

set(gca, 'FontName', 'Times New Roman');







