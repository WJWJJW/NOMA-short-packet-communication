% This script analysis blocklength and BLER by following step
% 1. solve common blocklength and alpha1 from e.q. 21 and 22 in ieee paper 
% under same theta value (figure 1)
% 2. check threshold : M1 thred and a1<1 thred
% 3. calculate optimal blocklength and alpha1 directly (figure 1)
% 4. plot blocklength v.s. BLER under optimal alpha1 (figure 2)
clc; clear variables; close all;
N = 1e6;
N1 = 256;
N2 = 256;


eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 20;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)


rho = pt/ no
% rho = db2pow(90);

% Set the user distance
d1 = 100;
d2 = 150;
eta = 4;

h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

lamda1 = mean(abs(h1).^2);
lamda2 = mean(abs(h2).^2);

delta = (300^-eta*eplsion2R)/(d1^-eta*eplsion1R);
contraint = (lamda1*eplsion1R)/(lamda2*eplsion2R);
theta = contraint*delta

alpha1 = 0.001:0.01:0.5;

% Analysis blcoklength and alpha1 under same theta value
for t = 1:length(alpha1)
    
    x1 = (1+theta*eplsion2R*lamda2*rho)/(1+theta*eplsion2R*lamda2*rho*alpha1(t));
    x2 = 1+alpha1(t)*rho*(lamda1*eplsion1R-theta*lamda2*eplsion2R);
    
    MD2(t) = N2/(log2 (x1));
    MD1(t) = N1/(log2 (x2));
end



if d2/d1 > (eplsion2R/eplsion1R)^0.25
    disp('constraint meet');
else
    disp('constraint not meet');  
end

if eplsion1R*lamda1 >= ((eplsion2R*lamda2*rho+4)*eplsion2R*lamda2)/(eplsion2R*lamda2*rho+2)
    disp('constaint 2 meet')
else 
    disp('constaint 2 not meet')
end


% optimal alpha1 and blocklength calculation
oopt_a1 = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
          4*(theta*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-theta*lamda2*eplsion2R))) ...
          / (2*theta*lamda2*eplsion2R*rho*(lamda1*eplsion1R-theta*lamda2*eplsion2R));
oopt_M = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(theta*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-theta*lamda2*eplsion2R)))...
                     / (2*theta*lamda2*eplsion2R)));

% Blocklength v.s. BLER under optimal alpha1                 
a = oopt_a1;
MM = 100:10:4000;

for m = 1:length(MM)
    BLER1(m) = (2^(N2/ MM(m))-1)/(lamda1*rho*(1-a-a*(2^(N2/ MM(m))-1)))+(2^(N1/MM(m))-1)/(lamda1*a*rho);
    BLER2(m) = (2^(N2/MM(m))-1)/(lamda2*rho*(1-a-a*(2^(N2/MM(m))-1)));
end


figure(1)
plot(alpha1, MD2, 'r', 'linewidth', 1.5);
hold on; grid on;
plot(alpha1, MD1, 'b', 'linewidth', 1.5);
plot(oopt_a1,oopt_M,'og','linewidth', 1.5);


axis([0 0.5 0 1000]);


xlabel('alpha1');
ylabel('M');
title('Required blocklength vs alpha1');
legend('M2','M1');

figure (2)
semilogy(MM, BLER2, 'r', 'linewidth', 1.5);
hold on; grid on;
semilogy(MM, BLER1, 'b', 'linewidth', 1.5);
xlabel('Blocklength');
ylabel('BLER');
title('BLER vs Blocklength');
legend('Far user','Near user');
xline(oopt_M,'-');
yline(1e-4,'-');
yline(1e-5,'-');

