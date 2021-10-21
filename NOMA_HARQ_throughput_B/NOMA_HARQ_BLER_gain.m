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
lambda1 = mean(abs(h1).^2);
lambda2 = mean(abs(h2).^2);

% Transmit power in dBm
Pt = 30;                
% Transmit power in linear scale
pt = (10^-3)*10^(Pt/10);

% AWGN noise
% BW = 10^6;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -80;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt/ no;

% # of information bits
N1 = 256;
N2 = 256;

% # of channel use (blocklength)
tx_times = 1:20;
M = 200;

% Power allocaitn coefficient
alphan = 0.1;
alphaf = 1 - alphan;

% BLER for NOMA-HARQ-IR
ir_eplsion_f = (2.^(N2./(tx_times.*M))-1)./(lambda2*rho*(alphaf-alphan.*(2.^(N1./(tx_times.*M))-1)));
ir_eplsion_nf = (2.^(N2./(tx_times.*M))-1)./(lambda1*rho*(alphaf-alphan.*(2.^(N1./(tx_times.*M))-1)));
ir_eplsion_nn = (2.^(N1./(tx_times.*M))-1)./(lambda1*rho*alphan);
ir_eplsion_n = ir_eplsion_nn + ir_eplsion_nf;

% BLER for NOMA-HARQ-CC
cc_eplsion_f = (2.^(N2./M)-1)./(lambda2*rho.*(tx_times.*alphaf-alphan.*(2.^(N1./M)-1)));
cc_eplsion_nf = (2.^(N2./M)-1)./(lambda1*rho.*(tx_times.*alphaf-alphan.*(2.^(N1./M)-1)));
cc_eplsion_nn = (2.^(N1./M)-1)./(lambda1*rho.*tx_times.*alphan);
cc_eplsion_n = cc_eplsion_nn + cc_eplsion_nf;

% Gain
Gain_f = ir_eplsion_f ./ cc_eplsion_f;
Gain_n = ir_eplsion_n ./ cc_eplsion_n;

figure (1);
plot(tx_times, Gain_f, 'b', 'linewidth', 1.5);
hold on; grid on;
plot(tx_times, Gain_n, 'r', 'linewidth', 1.5);

xlabel('Times of transmission');
ylabel('Gain');
legend('Far user','Near user');
set(gca, 'FontName', 'Times New Roman');

figure (2)
semilogy(tx_times, ir_eplsion_f, 'b', 'linewidth', 1.5);
hold on; grid on;
semilogy(tx_times, cc_eplsion_f, '-ob', 'linewidth', 1.5);
semilogy(tx_times, ir_eplsion_n, 'r', 'linewidth', 1.5);
semilogy(tx_times, cc_eplsion_n, '-or', 'linewidth', 1.5);

xlabel('Times of transmission');
ylabel('BLER');
legend('Far user w/ HARQ-IR','Far user w/ HARQ-cc','Near user w/ HARQ-IR','Near user w/ HARQ-CC');
set(gca, 'FontName', 'Times New Roman');


