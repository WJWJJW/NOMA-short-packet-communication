clc; clear variables; close all;
N = 1e6;
N1 = 256;
N2 = 256;


eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 0;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)


rho = pt/ no
% rho = db2pow(90);


d1 = 20;
d2 = 50;
eta = 4;

% lamda1 = 1/d1^eta;
% lamda2 = 1/d2^eta;

h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));

lamda1 = mean(abs(h1).^2);
lamda2 = mean(abs(h2).^2);

% t_delta1 = (lamda1*eplsion1R*(lamda2*eplsion2R*rho+2))/(2*lamda2*eplsion2R)-1
% t_delta2 = t_delta1/100

contraint = (lamda1*eplsion1R)/(lamda2*eplsion2R)
k = contraint/6;

alpha1 = 0.001:0.01:0.5;


for t = 1:length(alpha1)
    
    if d2/d1 > (eplsion2R/eplsion1R)^0.25
    x1 = (1+eplsion2R*lamda2*rho)/(1+eplsion2R*lamda2*rho*alpha1(t));
    x2 = 1+alpha1(t)*rho*(lamda1*eplsion1R-lamda2*eplsion2R);
    else
    x1 = (1+k*eplsion2R*lamda2*rho)/(1+k*eplsion2R*lamda2*rho*alpha1(t));
%     x2 = 1+(lamda1*eplsion1R*alpha1(t)*rho)/(t_delta2+1);
    x2 = 1+alpha1(t)*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R);
    end
    

    MD2(t) = N2/(log2 (x1));
    MD1(t) = N1/(log2 (x2));


%     syms M2
%     eqn2 = eplsion2R == (2^(N2/M2)-1)/(lamda2*rho*(1-alpha1(t)-alpha1(t)*(2^(N2/M2)-1)));
%     M2 = solve(eqn2,M2);
%     M21(t) = double (M2);
%     
%     syms M1
%     eqn1 = eplsion1R == (2^(N2/ M21(t))-1)/(lamda1*rho*(1-alpha1(t)-alpha1(t)*(2^(N2/ M21(t))-1)))+(2^(N1/M1)-1)/(lamda1*alpha1(t)*rho); 
%     eqn2 = eplsion1R == (2^(N2/ M1)-1)/(lamda1*rho*(1-alpha1(t)-alpha1(t)*(2^(N2/ M1)-1)))+(2^(N1/M1)-1)/(lamda1*alpha1(t)*rho); 
%     MA1 = solve(eqn1,M1);
%     MA11(t) = double (MA1 (1,1));
%     MB1 = solve(eqn2,M1);
%     MB11(t) = double (MB1 (1,1));



end



thred = (eplsion2R/eplsion1R)^0.25
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



% Solve alpha1 by Matlab
% syms otp_a1
% xo2 = (1+eplsion2R*lamda2*rho)/(1+eplsion2R*lamda2*rho*otp_a1);
% xo1 = 1+otp_a1*rho*(lamda1*eplsion1R-lamda2*eplsion2R);
% 
% eqn_opt = (xo2) == (xo1);
% opt_a1 = vpasolve(eqn_opt,otp_a1);


% optimal alpha1 and blocklength calculation
if d2/d1 > (eplsion2R/eplsion1R)^0.25
oopt_a1 = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
          4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R))) ...
          / (2*lamda2*eplsion2R*rho*(lamda1*eplsion1R-lamda2*eplsion2R));
oopt_M = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R)))...
                     / (2*lamda2*eplsion2R)));
else
% oopt_a1 = (-lamda1*eplsion1R - (t_delta2+1)*lamda2*eplsion2R + ...
%           sqrt((lamda1*eplsion1R + (t_delta2+1)*lamda2*eplsion2R)^2 ...
%           + 4*lamda1*eplsion1R*(lamda2*eplsion2R)^2*(t_delta2+1)*rho)) ...
%           / (2*lamda1*eplsion1R*lamda2*eplsion2R*rho)
% oopt_M = N1/log2(1+(lamda1*eplsion1R*oopt_a1*rho)/(t_delta2+1))

oopt_a1 = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
          4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R))) ...
          / (2*k*lamda2*eplsion2R*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R));
oopt_M = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R)))...
                     / (2*k*lamda2*eplsion2R)));
end


a = oopt_a1;
MM = 100:10:4000;

for m = 1:length(MM)
%     if d2/d1 > (eplsion2R/eplsion1R)^0.25
%     BLER1(m) = (2^(N2/ MM(m))-1)/(lamda1*rho*(1-a-a*(2^(N2/ MM(m))-1)))+(2^(N1/MM(m))-1)/(lamda1*a*rho);
%     else
%     BLER1(m) = (t_delta2+1)*(2^(N1/MM(m))-1)/(lamda1*a*rho);
%     end
    BLER1(m) = (2^(N2/ MM(m))-1)/(lamda1*rho*(1-a-a*(2^(N2/ MM(m))-1)))+(2^(N1/MM(m))-1)/(lamda1*a*rho);
    BLER2(m) = (2^(N2/MM(m))-1)/(lamda2*rho*(1-a-a*(2^(N2/MM(m))-1)));
end


figure(1)
plot(alpha1, MD2, 'r', 'linewidth', 1.5);
hold on; grid on;
plot(alpha1, MD1, 'b', 'linewidth', 1.5);
plot(oopt_a1,oopt_M,'og','linewidth', 1.5);

% plot(alpha1, M21,'or','linewidth', 1.5);
% plot(alpha1, MA11,'ob', 'linewidth', 1.5);
% plot(alpha1, MB11,'*g', 'linewidth', 1.5);
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







