% This script analysis relation between power and BLER
% fixed blocklength and power allocation coefficient

clc; clear variables; close all;
N1=80;
N2=80;
epsilon1R=10^-5;
epsilon2R=10^-4;


d1 = 100;
d2 = 170;

eta = 4;


h1 = sqrt(1/2*d1^-eta)*(randn(1,10000)+1i*randn(1,10000));
h2 = sqrt(1/2*d2^-eta)*(randn(1,10000)+1i*randn(1,10000));

lamda1 = mean(abs(h1).^2);
lamda2 = mean(abs(h2).^2);

% lamda1 = 1/(1+d1^eta);
% lamda2 = 1/(1+d2^eta);

BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

alpha1 = 0.499;
alpha2 = 1 - alpha1;
M = 100;

w1 = 2^(N1/M)-1;
w2 = 2^(N2/M)-1;

Xi_1 = sqrt(1/(2*pi*(2^(2*N1/M)-1)));
Xi_2 = sqrt(1/(2*pi*(2^(2*N2/M)-1)));

nu_1 = w1 - (1/(2*sqrt(M)*Xi_1));
tau_1 = w1 + (1/(2*sqrt(M)*Xi_1));
nu_2 = w2 - (1/(2*sqrt(M)*Xi_2));
tau_2 = w2 + (1/(2*sqrt(M)*Xi_2));

Pt = -10:30; %Transmitted power in dBm


p = length(Pt);
for u=1:p
%     rho = (10^-3)*10.^(RHO_db(u)/10);
%     rho = 10.^(RHO_db(u)/10);
    rho = (10^-3)*db2pow(Pt(u))/ no;


    % BLER from approximation Q + Riemann integral
    epsilon2(u) = 1 - exp(-w2/(lamda2*rho*(alpha2-alpha1*w2)));
    epsilon12(u) = 1 - exp(-w2/(lamda1*rho*(alpha2-alpha1*w2)));
    epsilon11(u) = 1 - exp(-w1/(lamda1*rho*alpha1));
    epsilon1(u) = epsilon12(u) + epsilon11(u);
    
    
    epsilon2_high_SNR(u) = w2/(lamda2*rho*(alpha2-alpha1*w2));
    epsilon12_high_SNR(u) = w2/(lamda1*rho*(alpha2-alpha1*w1));
    epsilon11_high_SNR(u) = w1/(lamda1*rho*alpha1);
    epsilon1_high_SNR(u) = epsilon12_high_SNR(u) + epsilon11_high_SNR(u);
    
    % BLER from approximation Q
    
    first_term_22 = ei(-alpha2/(lamda2*rho*alpha1*(alpha2-alpha1*nu_2)))-...
    ei(-alpha2/(lamda2*rho*alpha1*(alpha2-alpha1*tau_2)));
    
    first_term_12 = ei(-alpha2/(lamda1*rho*alpha1*(alpha2-alpha1*nu_2)))-...
    ei(-alpha2/(lamda1*rho*alpha1*(alpha2-alpha1*tau_2)));
    
    second_term_22 = exp(nu_2/(lamda2*rho*(alpha2-alpha1*nu_2)))*(alpha2-alpha1*nu_2) -...
    exp(nu_2/(lamda2*rho*(alpha2-alpha1*tau_2)))*(alpha2-alpha1*tau_2);
    
    second_term_12 = exp(nu_2/(lamda1*rho*(alpha2-alpha1*nu_2)))*(alpha2-alpha1*nu_2) - ...
    exp(nu_2/(lamda1*rho*(alpha2-alpha1*tau_2)))*(alpha2-alpha1*tau_2); 
    
    Epsilon2(u) = 1 - ((alpha2*Xi_2*sqrt(M)*exp(1/(lamda2*rho*alpha1)))/(lamda2*rho*(alpha1^2)))* first_term_22...
    - ((Xi_2*sqrt(M))/alpha1)*second_term_22;
    
    Epsilon12(u) = 1 - ((alpha2*Xi_2*sqrt(M)*exp(1/(lamda1*rho*alpha1)))/(lamda1*rho*(alpha1^2)))* first_term_12...
    - ((Xi_2*sqrt(M))/alpha1)*second_term_12;
    
    
    Epsilon11(u) = 1 + Xi_1*sqrt(M)*lamda1*alpha1*rho*(exp(-tau_1/(lamda1*alpha1*rho))-exp(-nu_1/(lamda1*alpha1*rho)));
    Epsilon1(u) = Epsilon12(u) + Epsilon11(u);

    % BLER from approximation Q by MATLAB
    fun1 = @(x) 1-exp(x./(-lamda2*rho*(alpha2-alpha1*x)));
    fun2 = @(x) 1-exp(x./(-lamda1*rho*(alpha2-alpha1*x)));
    fun3 = @(x) 1-exp(x./(-lamda1*rho*alpha1));

    MAT_Epsilon2(u) = Xi_2*sqrt(M)*integral(fun1,nu_2,tau_2);
    MAT_Epsilon12(u) = Xi_2*sqrt(M)*integral(fun2,nu_2,tau_2);
    MAT_Epsilon11(u) = Xi_1*sqrt(M)*integral(fun3,nu_1,tau_1);
    MAT_Epsilon1(u) = MAT_Epsilon12(u) + MAT_Epsilon11(u);
    
end

semilogy(Pt, epsilon2,'or'); hold on; grid on;
semilogy(Pt, epsilon1, 'ob');
semilogy(Pt, epsilon12,'og');
semilogy(Pt, epsilon11, 'ok');

% semilogy(Pt, Epsilon2,'+r', 'linewidth',2);
% semilogy(Pt, Epsilon1, '+b', 'linewidth',2);
% semilogy(Pt, Epsilon12,'+g', 'linewidth',2);
% semilogy(Pt, Epsilon11, '+k', 'linewidth',2);

semilogy(Pt, MAT_Epsilon2,'r', 'linewidth',1.5);
semilogy(Pt, MAT_Epsilon1, 'b','linewidth',1.5);
semilogy(Pt, MAT_Epsilon12,'g', 'linewidth',1.5);
semilogy(Pt, MAT_Epsilon11, 'k','linewidth',1.5);

semilogy(Pt, epsilon2_high_SNR,'*r');
semilogy(Pt, epsilon1_high_SNR, '*b');
semilogy(Pt, epsilon12_high_SNR,'*g');
semilogy(Pt, epsilon11_high_SNR, '*k');

title('BLER vs Transmit Power');
xlabel('Transmit power (P in dBm)');
ylabel('BLER');
% legend('epsilon2','epsilon1','epsilon12','epsilon11','Epsilon2','Epsilon1','Epsilon12','Epsilon11',...
%        'MAT Epsilon2','MAT Epsilon1','MAT Epsilon12','MAT Epsilon11',...
%        'epsilon2 high SNR', 'epsilon1 high SNR', 'epsilon12 high SNR', 'epsilon11 high SNR' );

legend('Q+R-appor.2','Q+R-appor.1','Q+R-appor.12','Q+R-appor.11',...
       'Q-appor.2','Q-appor.1','Q-appor.12','Q-appor.11',...
       'epsilon2 high SNR', 'epsilon1 high SNR', 'epsilon12 high SNR', 'epsilon11 high SNR' );

