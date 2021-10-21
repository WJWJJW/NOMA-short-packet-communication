clc; clear variables; close all;

N = 1e6; % number of channel tap
NN = 256; % number of information bit
N1 = NN;
N2 = NN;

% Target BLER
eplsion1R = 10^-5;
eplsion2R = 10^-4;


Pt = 30;               %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

% AWGN
% BW = 10^7;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt / no;

% Set the user distance
D1 = [100 150];
d2 = 200;
eta = 4;

delta = 0.2:0.001:0.7;
opt_delta = zeros(1,length(D1));
Opt_M_delta = zeros(1,length(D1));


dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);

Opt_a1 = zeros(length(Pt),length(delta)); 
Opt_M = zeros(length(Pt),length(delta));

for u=1:length(D1)
    d1 = D1(u);
    % Rayleigh fading channel
    h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
    h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

    % Channel mean
    lamda1 = mean(abs(h1).^2);
    lamda2 = mean(abs(h2).^2);
    contraint = (lamda1*eplsion1R)/(lamda2*eplsion2R);
    for dd = 1:length(delta)
        theta = contraint*delta(dd);
        check = (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - ...
            (lamda2*eplsion2R*rho + 4)*lamda2*eplsion2R > 0;
        switch (check)
            case 0
                eplsion2R_m = theta*eplsion2R;
            case 1
                eplsion2R_m = eplsion2R;
        end
        Opt_a1(u,dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
        4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m))) ...
        / (2*lamda2*eplsion2R_m*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m));
    
        Opt_M(u,dd) = N1/log2(1+((-lamda1*eplsion1R + ...
        sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)))...
        / (2*lamda2*eplsion2R_m))); 
    end
    % Find optimal delta 
    [opt_delta(u)] = delta_finder (lamda1*rho*eplsion1R);
    % Optimal Blocklength
    theta = contraint*opt_delta(u);
        check = (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - ...
            (lamda2*eplsion2R*rho + 4)*lamda2*eplsion2R > 0;
        switch (check)
            case 0
                eplsion2R_m = theta*eplsion2R;
            case 1
                eplsion2R_m = eplsion2R;
        end
    
        Opt_M_delta(u) = N1/log2(1+((-lamda1*eplsion1R + ...
        sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)))...
        / (2*lamda2*eplsion2R_m))); 
    
end




figure (1)
yyaxis left
plot(delta, Opt_a1(1,:),'b','linewidth',1.5);
hold on; grid on;
plot(delta, Opt_a1(2,:), '--b','linewidth',1.5);

yline(0.5,'-','HandleVisibility','off');
ylabel('\alpha_1');

yyaxis right
plot(delta, Opt_M(1,:), 'Color',[1 0.5 0],'linewidth',1.5);
hold on; grid on;
plot(delta, Opt_M(2,:),'--', 'Color',[1 0.5 0],'linewidth',1.5);

plot(opt_delta(1), Opt_M_delta(1),'og','linewidth',1.5);
plot(opt_delta(2), Opt_M_delta(2),'om','linewidth',1.5);

xlabel('\delta')
ylabel('Blocklength (Channel uses)');
legend('\alpha_i^{opt} for d_i=100m','\alpha_i^{opt} for d_i=150m',...
        'M_{opt}^{1 pair} for d_i=100m', 'M_{opt}^{1 pair} for d_i=150m');
    
    
ar = annotation('arrow');
x = [0.3 0.5];
y = [0.6 0.5];
an = annotation('textarrow',x,y,'String',...
    '\delta_{opt} by (41)');
an.FontName = 'Times New Roman';

set(gca, 'FontName', 'Times New Roman'); 

