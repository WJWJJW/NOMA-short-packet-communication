clear all;close all;clc;

N = 10000;
eta = 4;
% User distance
d1 = 400;
d2 = 700;

% Rayleigh fading channel
h1 = sqrt(1/2*d1^-eta)*(randn(N,1)+1i*randn(N,1));
h2 = sqrt(1/2*d2^-eta)*(randn(N,1)+1i*randn(N,1));

% Transmit power in dBm
Pt = 0:2:50;                
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);

% AWGN noise
BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)



%%%%% Fixed transmitted power %%%%%

alpha1 = 0.01:0.01:0.5;
alpha2 = 1- alpha1;

% Received SINR
gamma2(1,:) = mean((alpha2*pt(length(Pt)).*abs(h2).^2)./(alpha1*pt(length(Pt)).*abs(h2).^2+no));
gamma12(1,:) = mean((alpha2*pt(length(Pt)).*abs(h1).^2)./(alpha1*pt(length(Pt)).*abs(h1).^2+no));

gamma2(2,:) = mean((alpha2*pt(11).*abs(h2).^2)./(alpha1*pt(11).*abs(h2).^2+no));
gamma12(2,:) = mean((alpha2*pt(11).*abs(h1).^2)./(alpha1*pt(11).*abs(h1).^2+no));

gamma2(3,:) = mean((alpha2*pt(1).*abs(h2).^2)./(alpha1*pt(1).*abs(h2).^2+no));
gamma12(3,:) = mean((alpha2*pt(1).*abs(h1).^2)./(alpha1*pt(1).*abs(h1).^2+no));

figure (1)
subplot(3,1,1);
plot(alpha1, alpha2 - alpha1.*gamma2(1,:), 'linewidth', 1.5);
hold on; grid on;
plot(alpha1, alpha2 - alpha1.*gamma12(1,:), 'linewidth', 1.5);
title('$$\alpha_2$$ - $$\alpha_1$$ $$\times$$ $$\gamma $$ vs $$\alpha_1$$', 'Interpreter', 'latex');
xlabel('$\alpha_1$', 'Interpreter', 'latex');
legend('\gamma_{2} w/ P=50dBm','\gamma_{12} w/ P=50dBm');
subplot(3,1,2);
plot(alpha1, alpha2 - alpha1.*gamma2(2,:), 'linewidth', 1.5);
hold on; grid on;
plot(alpha1, alpha2 - alpha1.*gamma12(2,:), 'linewidth', 1.5);
title('$$\alpha_2$$ - $$\alpha_1$$ $$\times$$ $$\gamma $$ vs $$\alpha_1$$', 'Interpreter', 'latex');
xlabel('$\alpha_1$', 'Interpreter', 'latex');
legend('\gamma_{2} w/ P=20dBm','\gamma_{12} w/ P=20dBm');
subplot(3,1,3);
plot(alpha1, alpha2 - alpha1.*gamma2(3,:), 'linewidth', 1.5);
hold on; grid on;
plot(alpha1, alpha2 - alpha1.*gamma12(3,:), 'linewidth', 1.5);
title('$$\alpha_2$$ - $$\alpha_1$$ $$\times$$ $$\gamma $$ vs $$\alpha_1$$', 'Interpreter', 'latex');
xlabel('$\alpha_1$', 'Interpreter', 'latex');
legend('\gamma_{2} w/ P=0dBm','\gamma_{12} w/ P=0dBm');






%%%%% Fixed power allocation coefficient %%%%%
a1 = 0.4;
a2 = 1- a1;

for u = 1:length(Pt)
    % Received SINR
    Gamma2(u) = mean((a2*pt(u).*abs(h2).^2)./(a1*pt(u).*abs(h2).^2+no));
    check2(u) = a2 - a1.*Gamma2(u);
    Gamma12(u) = mean((a2*pt(u).*abs(h1).^2)./(a1*pt(u).*abs(h1).^2+no));
    check12(u) = a2 - a1.*Gamma12(u);
end

figure (2)
plot(Pt, check2, 'linewidth', 1.5);
hold on; grid on;
plot(Pt, check12, 'linewidth', 1.5);

title('$$\alpha_2$$ - $$\alpha_1$$ $$\times$$ $$\gamma $$ vs Transmit Power', 'Interpreter', 'latex');
xlabel('Transmit power (in dBm)');
legend('\gamma_{2}','\gamma_{12}');

%%%%% Fixed power allocation coefficient and transmitted power Part 1 %%%%%

p = (10^-3)*10.^(40/10);

% User distance
d1 = 400;
d_diff = 50:10:600;

for rr = 1:length(d_diff)
    d2 = d1+d_diff(rr);
    % Rayleigh fading channel
    hd1 = sqrt(1/2*d1^-eta)*(randn(N,1)+1i*randn(N,1));
    hd2 = sqrt(1/2*d2^-eta)*(randn(N,1)+1i*randn(N,1));
    
    % Received SINR
    ggamma2(rr) = mean((a2*p.*abs(hd2).^2)./(a1*p.*abs(hd2).^2+no));
    ccheck2(rr) = a2 - a1.*ggamma2(rr);
    ggamma12(rr) = mean((a2*p.*abs(hd1).^2)./(a1*p.*abs(hd1).^2+no));
    ccheck12(rr) = a2 - a1.*ggamma12(rr);
end

figure (3)
plot(d_diff, ccheck2, 'linewidth', 1.5);
hold on; grid on;
plot(d_diff, ccheck12, 'linewidth', 1.5);

title('$$\alpha_2$$ - $$\alpha_1$$ $$\times$$ $$\gamma $$ vs Distance difference', 'Interpreter', 'latex');
xlabel('Distance difference (m)');
legend('\gamma_{2}','\gamma_{12}');

%%%%% Fixed power allocation coefficient and transmitted power Part 2 %%%%%
D1 = [100, 200, 300, 400, 500, 600, 700];
D2 = [200, 400, 600, 800, 1000, 1200, 1400];

for rrr = 1:length(D1)
    
    % Rayleigh fading channel
    hD1 = sqrt(1/2*D1(rrr)^-eta)*(randn(N,1)+1i*randn(N,1));
    hD2 = sqrt(1/2*D2(rrr)^-eta)*(randn(N,1)+1i*randn(N,1));
    
    % Received SINR
    gggamma2(rrr) = mean((a2*p.*abs(hD2).^2)./(a1*p.*abs(hD2).^2+no));
    cccheck2(rrr) = a2 - a1.*gggamma2(rrr);
    gggamma12(rrr) = mean((a2*p.*abs(hD1).^2)./(a1*p.*abs(hD1).^2+no));
    cccheck12(rrr) = a2 - a1.*gggamma12(rrr);
end

rrr = 1:length(D1);
figure (4)
semilogy(rrr, cccheck2, 'linewidth', 1.5);
hold on; grid on;
semilogy(rrr, cccheck12, 'linewidth', 1.5);

title('$$\alpha_2$$ - $$\alpha_1$$ $$\times$$ $$\gamma $$ vs Distance set', 'Interpreter', 'latex');
xlabel('Distance Set');
legend('\gamma_{2}','\gamma_{12}');





