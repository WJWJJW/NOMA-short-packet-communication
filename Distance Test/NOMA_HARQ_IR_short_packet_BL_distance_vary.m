clc; clear variables; close all;
N = 1e6;
NN = 80;
N1 = NN;
N2 = NN;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 20;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt/ no;
RHO = pow2db(rho);

eta = 4;
d1 = 10:1:150;
d2 = 10:1:400;

opt_a1 = zeros(length(d1),length(d2));
opt_M = zeros(length(d1),length(d2));
h = randn(1,N)+1i*randn(1,N);

for d=1:length(d1)
    for dd=1:length(d2)
        h1 = sqrt(1/2*d1(d)^-eta).*h;
        h2 = sqrt(1/2*d2(dd)^-eta).*h;

        lamda1 = mean(abs(h1).^2);
        lamda2 = mean(abs(h2).^2);
        contraint(d,dd) = (lamda1*eplsion1R)/(lamda2*eplsion2R);
        k = contraint(d,dd)/3;


        % optimal alpha1 and blocklength calculation
        if d2(dd)/d1(d) > (eplsion2R/eplsion1R)^0.25 && (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - (lamda2*eplsion2R.*rho + 4)*lamda2*eplsion2R > 0
            opt_a1(d,dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
                      4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R))) ...
                      / (2*lamda2*eplsion2R*rho*(lamda1*eplsion1R-lamda2*eplsion2R));
            opt_M(d,dd) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-lamda2*eplsion2R)))...
                                 / (2*lamda2*eplsion2R)));
        elseif d2(dd)/d1(d) > (eplsion2R/eplsion1R)^0.25 && (lamda2*eplsion2R*rho + 2)*lamda1*eplsion1R - (lamda2*eplsion2R.*rho + 4)*lamda2*eplsion2R < 0
            opt_a1(d,dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
                      4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R))) ...
                      / (2*k*lamda2*eplsion2R*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R));
            opt_M(d,dd) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R)))...
                                 / (2*k*lamda2*eplsion2R)));
        else
            opt_a1(d,dd) = (-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+...
                      4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R))) ...
                      / (2*k*lamda2*eplsion2R*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R));
            opt_M(d,dd) = N1/log2(1+((-lamda1*eplsion1R + sqrt((lamda1*eplsion1R)^2+4*(k*lamda2*eplsion2R)^2*rho*(lamda1*eplsion1R-k*lamda2*eplsion2R)))...
                                 / (2*k*lamda2*eplsion2R)));
        end
    end

end

dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);

figure (1)
yyaxis left
plot(d2, real(opt_a1(1,:)),'-b','linewidth', 1.5);
hold on; grid on;
plot(d2, real(opt_a1(10,:)),'-g','linewidth', 1.5);
plot(d2, real(opt_a1(50,:)),'-m','linewidth', 1.5);
plot(d2, real(opt_a1(100,:)),'-r','linewidth', 1.5);
plot(d2, real(opt_a1(140,:)),'-c','linewidth', 1.5);


yline(0.5,'-');

ylabel('\alpha_1');


yyaxis right
plot(d2, real(opt_M(1,:)),'ob');
hold on; grid on;
plot(d2, real(opt_M(10,:)),'og');
plot(d2, real(opt_M(50,:)),'om');
plot(d2, real(opt_M(100,:)),'or');
plot(d2, real(opt_M(140,:)),'oc');
ylabel('blocklength');

xline(d1(1)*dis_thred1,'-','Threshold1a');
xline(d1(10)*dis_thred1,'-','Threshold1b');
xline(d1(50)*dis_thred1,'-','Threshold1c');
xline(d1(100)*dis_thred1,'-','Threshold1d');
xline(d1(140)*dis_thred1,'-','Threshold1e');


xlabel('Distance');

figure (2)
yyaxis left
plot(d2, real(opt_a1(1,:)),'-b','linewidth', 1.5);
hold on; grid on;




yline(0.5,'-');

ylabel('\alpha_1');


yyaxis right
plot(d2, real(opt_M(1,:)),'ob');
hold on; grid on;



ylabel('blocklength');

xline(d1(1)*dis_thred1,'-','Threshold1a');




figure (3)
yyaxis left
plot(d2, real(opt_a1(50,:)),'-b','linewidth', 1.5);
hold on; grid on;
plot(d2, real(opt_a1(70,:)),'-g','linewidth', 1.5);
plot(d2, real(opt_a1(80,:)),'-m','linewidth', 1.5);



yline(0.5,'-');

ylabel('\alpha_1');


yyaxis right
plot(d2, real(opt_M(50,:)),'ob');
hold on; grid on;
plot(d2, real(opt_M(70,:)),'og');
plot(d2, real(opt_M(80,:)),'om');


ylabel('blocklength');

xline(d1(50)*dis_thred1,'-','Threshold1a');
xline(d1(70)*dis_thred1,'-','Threshold1b');
xline(d1(80)*dis_thred1,'-','Threshold1c');




figure (4)

yyaxis left
plot(d2, real(opt_a1(100,:)),'-b','linewidth', 1.5);
hold on; grid on;
plot(d2, real(opt_a1(120,:)),'-g','linewidth', 1.5);
plot(d2, real(opt_a1(140,:)),'-m','linewidth', 1.5);


yline(0.5,'-');

ylabel('\alpha_1');


yyaxis right
plot(d2, real(opt_M(100,:)),'ob');
hold on; grid on;
plot(d2, real(opt_M(120,:)),'og');
plot(d2, real(opt_M(140,:)),'om');

ylabel('blocklength');

xline(d1(100)*dis_thred1,'-','Threshold1a');
xline(d1(120)*dis_thred1,'-','Threshold1b');
xline(d1(140)*dis_thred1,'-','Threshold1c');




