clc; clear variables; close all;
N1=80;
N2=80;
epsilon1R=10^-5;
epsilon2R=10^-4;

d1 = 400;
eta = 4;


BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

alpha1 = 0.1;
alpha2 = 1 - alpha1;
M = 100;

w1 = 2^(N1/M)-1;
w2 = 2^(N2/M)-1;

Xi_1 = sqrt(1/(2*pi*(2^(2*N1/M)-1)));
Xi_2 = sqrt(1/(2*pi*(2^(2*N2/M)-1)));

eta_1 = w1 - (1/(2*sqrt(M)*Xi_1));
tau_1 = w1 + (1/(2*sqrt(M)*Xi_1));
eta_2 = w2 - (1/(2*sqrt(M)*Xi_2));
tau_2 = w2 + (1/(2*sqrt(M)*Xi_2));

Pt = 50; %Transmitted power in dBm
rho = (10^-3)*db2pow(Pt)/ no;

d_diff = 50:50:450;
for u=1:length(d_diff)
    d2 = d1+d_diff(u);
    h1 = sqrt(1/2*d1^-eta)*(randn(1,10000)+1i*randn(1,10000));
    h2 = sqrt(1/2*d2^-eta)*(randn(1,10000)+1i*randn(1,10000));

    lamda1 = mean(abs(h1).^2);
    lamda2 = mean(abs(h2).^2);
    
    % BLER from approximation Q + Riemann integral
    epsilon2(u) = 1 - exp(-w2/(lamda2*rho*(alpha2-alpha1*w2)));
    epsilon12(u) = 1 - exp(-w2/(lamda1*rho*(alpha2-alpha1*w2)));
    epsilon11(u) = 1 - exp(-w1/(lamda1*rho*alpha1));
    epsilon1(u) = epsilon12(u) + epsilon11(u);
    
    % BLER from approximation Q + Riemann integral + high SNR approximation
    epsilon2_high_SNR(u) = w2/(lamda2*rho*(alpha2-alpha1*w2));
    epsilon12_high_SNR(u) = w2/(lamda1*rho*(alpha2-alpha1*w1));
    epsilon11_high_SNR(u) = w1/(lamda1*rho*alpha1);
    epsilon1_high_SNR(u) = epsilon12_high_SNR(u) + epsilon11_high_SNR(u);
    
    % BLER from approximation Q by MATLAB
    fun1 = @(x) 1-exp(x./(-lamda2*rho*(alpha2-alpha1*x)));
    fun2 = @(x) 1-exp(x./(-lamda1*rho*(alpha2-alpha1*x)));
    fun3 = @(x) 1-exp(x./(-lamda1*rho*alpha1));

    MAT_Epsilon2(u) = Xi_2*sqrt(M)*integral(fun1,eta_2,tau_2);
    MAT_Epsilon12(u) = Xi_2*sqrt(M)*integral(fun2,eta_2,tau_2);
    MAT_Epsilon11(u) = Xi_1*sqrt(M)*integral(fun3,eta_1,tau_1);
    MAT_Epsilon1(u) = MAT_Epsilon12(u) + MAT_Epsilon11(u);
end




semilogy(d_diff, epsilon2,'or'); hold on; grid on;
semilogy(d_diff, epsilon1, 'ob');
semilogy(d_diff, epsilon12,'og');
semilogy(d_diff, epsilon11, 'ok');

semilogy(d_diff, MAT_Epsilon2,'r', 'linewidth',1.5);
semilogy(d_diff, MAT_Epsilon1, 'b','linewidth',1.5);
semilogy(d_diff, MAT_Epsilon12,'g', 'linewidth',1.5);
semilogy(d_diff, MAT_Epsilon11, 'k','linewidth',1.5);

semilogy(d_diff, epsilon2_high_SNR,'*r');
semilogy(d_diff, epsilon1_high_SNR, '*b');
semilogy(d_diff, epsilon12_high_SNR,'*g');
semilogy(d_diff, epsilon11_high_SNR, '*k');

title('BLER vs Difference of distance');
xlabel('Difference of distance (in meter)');
ylabel('BLER');


legend('epsilon2','epsilon1','epsilon12','epsilon11',...
       'MAT Epsilon2','MAT Epsilon1','MAT Epsilon12','MAT Epsilon11',...
       'epsilon2 high SNR', 'epsilon1 high SNR', 'epsilon12 high SNR', 'epsilon11 high SNR' );

