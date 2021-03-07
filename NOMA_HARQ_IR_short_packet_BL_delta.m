clc; clear variables; close all;

N = 1e6; % number of channel tap
NNN = 1000; % number of Monte Carlo
K = 5;  % number of cluster (number of user  = 2K)
NN = 256; % number of information bit
N1 = NN;
N2 = NN;

eplsion1R = 10^-5;
eplsion2R = 10^-4;


Pt = 20:1:30;               %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

% AWGN
% BW = 10^7;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

% Set the user distance
d1 = 100;
d2 = 170;
eta = 4;

% Rayleigh fading channel
h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

% Channel mean
lamda1 = mean(abs(h1).^2);
lamda2 = mean(abs(h2).^2);

delta = [1/4 1/8];

dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);

Opt_M = zeros(length(Pt),length(Pt)+1);

for u=1:length(Pt)
    rho = pt(u) / no;
    % Find optimal delta 
    [opt_delta] = delta_finder (lamda1*rho*eplsion1R);
    contraint = (lamda1*eplsion1R)/(lamda2*eplsion2R);
    for dd = 1: length(delta)
        theta = contraint*delta(dd);
        check = (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - ...
            (lamda2*eplsion2R*rho + 4)*lamda2*eplsion2R > 0;
        switch (check)
            case 0
                eplsion2R_m = theta*eplsion2R;
            case 1
                eplsion2R_m = eplsion2R;
        end
        Opt_M(u,dd) = N1/log2(1+((-lamda1*eplsion1R + ...
        sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)))...
        / (2*lamda2*eplsion2R_m))); 
    end
    theta = contraint*opt_delta;
    check = (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - ...
        (lamda2*eplsion2R*rho + 4)*lamda2*eplsion2R > 0;
    switch (check)
        case 0
            eplsion2R_m = theta*eplsion2R;
        case 1
            eplsion2R_m = eplsion2R;
    end
    Opt_M(u,3) = N1/log2(1+((-lamda1*eplsion1R + ...
    sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R_m)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R_m)))...
    / (2*lamda2*eplsion2R_m))); 
end



figure (1)

plot(Pt, Opt_M(:,1),'b');
hold on; grid on;
plot(Pt, Opt_M(:,2), 'Color',[1 0.5 0]);
plot(Pt, Opt_M(:,3), 'r');


xlabel('Transmitted Power (dBm)')
ylabel('Blocklength (Channel uses)');
legend('\delta = 1/4','\delta = 1/8', '\delta by (41)');
set(gca, 'FontName', 'Times New Roman'); 
