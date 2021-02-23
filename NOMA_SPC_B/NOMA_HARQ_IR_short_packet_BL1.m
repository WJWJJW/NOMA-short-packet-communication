% This script analysis blocklength and BLER by following step
% 1. solve common blocklength and alpha1 from e.q. 21 and 22 in ieee paper 
% under diferent condition (theta value: 1 or other) (figure 1)
% 2. check threshold : M1 thred and a1<1 thred
% 3. calculate optimal blocklength and alpha1 directly (figure 1)
% 4. plot blocklength v.s. BLER under optimal alpha1 (figure 2)
clc; clear variables; close all;
% # of channel tap
N = 1e7;
% # of infomation bit
N1 = 256;
N2 = 256;

% Target BLER
eplsion1R = 10^-4;
eplsion2R = 10^-4;

% Transmit power in dBm
Pt = 30;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

% AWGN noise
% BW = 10^6;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

% Tx SNR
rho = pt/ no;

% Set the user distance
d1 = 100;
d2 = 120;
eta = 4;

% Rayleigh fading channel
h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

% Channel mean
lamda1 = mean(abs(h1).^2);
lamda2 = mean(abs(h2).^2);

% Condition check
if d2/d1 > (eplsion2R/eplsion1R)^0.25
    disp('Condition 1 meet');
else
    disp('Condition 1 not meet');  
end

if eplsion1R*lamda1 >= ((eplsion2R*lamda2*rho+4)*eplsion2R*lamda2)/(eplsion2R*lamda2*rho+2)
    disp('Condition 2 meet')
else 
    disp('Condition 2 not meet')
end

% Find optimal delta according to near user
[delta] = delta_finder (lamda1*rho*eplsion1R);
theta = delta*(lamda1*eplsion1R)/(lamda2*eplsion2R);

% Adjust target BLER according to condition 2
check = (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - (lamda2*eplsion2R.*rho + 4)*lamda2*eplsion2R > 0;
switch (check)
    case 0
        eplsion2R_m = theta*eplsion2R;
    case 1
        eplsion2R_m = eplsion2R;
end

alpha1 = 0.001:0.01:0.5;

% Analysis blcoklength and alpha1 under different theta value
for t = 1:length(alpha1)
    x1 = (1+eplsion2R_m*lamda2*rho)/(1+eplsion2R_m*lamda2*rho*alpha1(t));
    x2 = 1+alpha1(t)*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m);

    MD2(t) = N2/(log2 (x1));
    MD1(t) = N1/(log2 (x2));
end

% Optimal alpha1 and blocklength calculation
oopt_a1 = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
          4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m))) ...
          / (2*lamda2*eplsion2R_m*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m))
oopt_M = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)))...
                     / (2*lamda2*eplsion2R_m)))


% Blocklength v.s. BLER under optimal alpha1
a = oopt_a1;
MM = 150:10:3000;

for m = 1:length(MM)
    BLER1(m) = (2^(N2/ MM(m))-1)/(lamda1*rho*(1-a-a*(2^(N2/ MM(m))-1)))+(2^(N1/ MM(m))-1)/(lamda1*a*rho);
    BLER2(m) = (2^(N2/ MM(m))-1)/(lamda2*rho*(1-a-a*(2^(N2/ MM(m))-1)));
end


figure(1)
plot(alpha1, MD2, 'r', 'linewidth', 1.5);
hold on; grid on;
plot(alpha1, MD1, 'b', 'linewidth', 1.5);
plot(oopt_a1,oopt_M,'og','linewidth', 1.5);

axis([0 0.5 100 1500]);

xlabel('\alpha_1');
ylabel('Blocklength (Channel uses)');
% title('Required blocklength vs \alpha_1');
legend('Far user','Near user');

x = [0.3 0.5];
y = [0.6 0.5];
an = annotation('textarrow',x,y,'String',...
    '\alpha_1^{opt} and M_{opt}^{1 pair} by using (35) and (36)');
an.FontName = 'Times New Roman';

set(gca, 'FontName', 'Times New Roman');

figure (2)
semilogy(MM, BLER2, 'r', 'linewidth', 1.5);
hold on; grid on;
semilogy(MM, BLER1, 'b', 'linewidth', 1.5);
xlabel('Blocklength (Channel uses)');
ylabel('BLER');
% title('BLER vs Blocklength');
legend('Far user','Near user');
xx = xline(oopt_M,'-','Optimal Required Blocklength','HandleVisibility','off');
xx.FontName = 'Times New Roman';
yline(1e-4,'-','HandleVisibility','off');
yline(1e-5,'-','HandleVisibility','off');

set(gca, 'FontName', 'Times New Roman');

