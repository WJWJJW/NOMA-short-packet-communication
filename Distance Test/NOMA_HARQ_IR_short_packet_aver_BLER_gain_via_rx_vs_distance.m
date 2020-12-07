clc; clear variables; close all;
tic
%% Basic setting
N = 1e6;
eta = 4;
% User distance
d1 = 400;
d2 = 500:100:1500;

% d1 = 100:100:700;
% d2 = 1000;


% AWGN noise
BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

% Parameter setting
M1 = [100 200 300];
M2 = [100 200 300];
N1 = 80;
N2 = 80;

w1 = 2.^(N1./M1)-1;
w2 = 2.^(N2./M2)-1;

% Transmit power in dBm
Pt = 40;                
% Pt = -113:2:-74;
% Transmit power in linear scale
pt = (10^-3)*10^(Pt/10);

rho = (10^-3)*db2pow(Pt)/ no;
% Power allocaitn coefficient
alpha1 = 0.3;
alpha2 = 1 - alpha1;
%%
for mm = 1:length(M1)
    for dd = 1:length(d2)
        % Rayleigh fading channel
        h1 = sqrt(1/2*d1^-eta)*(randn(1,N)+1i*randn(1,N));
        h2 = sqrt(1/2*d2(dd)^-eta)*(randn(1,N)+1i*randn(1,N));
%         h1 = sqrt(1/2*d1(dd)^-eta)*(randn(1,N)+1i*randn(1,N));
%         h2 = sqrt(1/2*d2^-eta)*(randn(1,N)+1i*randn(1,N));
        
        lambda1 = mean(abs(h1).^2);
        lambda2 = mean(abs(h2).^2);
        
        %% Q-function w/ NO approximation
        % Received SINR
        gamma2(dd,:) = (alpha2*pt.*abs(h2).^2)./(alpha1*pt.*abs(h2).^2+no);
        gamma12(dd,:) = (alpha2*pt.*abs(h1).^2)./(alpha1*pt.*abs(h1).^2+no);
        gamma11(dd,:) = (alpha1*pt.*abs(h1).^2)./no;

        % Q-function argument
        x_2a(mm,dd,:) = log2 (1+gamma2(dd,:)).*M2(mm) - N2;
        x_2b(mm,dd,:) = log2 (exp(1))*sqrt((1-1./(1+gamma2(dd,:)).^2) .*M2(mm));
        ratio2(mm,dd,:) = x_2a(mm,dd,:)./x_2b(mm,dd,:);

        x_12a(mm,dd,:) = log2 (1+gamma12(dd,:)).*M1(mm) - N2;
        x_12b(mm,dd,:) = log2 (exp(1))*sqrt((1-1./(1+gamma12(dd,:)).^2) .*M1(mm));
        ratio12(mm,dd,:) = x_12a(mm,dd,:)./x_12b(mm,dd,:);

        x_11a(mm,dd,:) = log2 (1+gamma11(dd,:)).*M1(mm) - N1;
        x_11b(mm,dd,:) = log2 (exp(1))*sqrt((1-1./(1+gamma11(dd,:)).^2) .*M1(mm));
        ratio11(mm,dd,:) = x_11a(mm,dd,:)./x_11b(mm,dd,:);

        % Q-function
        q2(mm,dd,:) = qfunc(ratio2(mm,dd,:));
        q12(mm,dd,:) = qfunc(ratio12(mm,dd,:)); 
        q11(mm,dd,:) = qfunc(ratio11(mm,dd,:));
    end
end

% BLER1_gain(1,:) = mean(q12(1,:,:),3)+(1-mean(q12(1,:,:),3)).*mean(q11(1,:,:),3)-...
%                   (mean(q12(2,:,:),3)+(1-mean(q12(2,:,:),3)).*mean(q11(2,:,:),3));
% BLER1_gain(2,:) = mean(q12(2,:,:),3)+(1-mean(q12(2,:,:),3)).*mean(q11(2,:,:),3)-...
%                   (mean(q12(3,:,:),3)+(1-mean(q12(3,:,:),3)).*mean(q11(3,:,:),3));
% BLER1_gain(3,:) = mean(q12(1,:,:),3)+(1-mean(q12(1,:,:),3)).*mean(q11(1,:,:),3)-...
%                   (mean(q12(3,:,:),3)+(1-mean(q12(3,:,:),3)).*mean(q11(3,:,:),3));
%               
% BLER2_gain(1,:) = mean(q2(1,:,:),3) - mean(q2(2,:,:),3);
% BLER2_gain(2,:) = mean(q2(2,:,:),3) - mean(q2(3,:,:),3);
% BLER2_gain(3,:) = mean(q2(1,:,:),3) - mean(q2(3,:,:),3);


BLER1_gain(1,:) = pow2db((mean(q12(2,:,:),3)+(1-mean(q12(2,:,:),3)).*mean(q11(2,:,:),3) )./...
                  (mean(q12(1,:,:),3)+(1-mean(q12(1,:,:),3)).*mean(q11(1,:,:),3)));
BLER1_gain(2,:) = pow2db((mean(q12(3,:,:),3)+(1-mean(q12(3,:,:),3)).*mean(q11(3,:,:),3)) ./...
                  (mean(q12(2,:,:),3)+(1-mean(q12(2,:,:),3)).*mean(q11(2,:,:),3)));
BLER1_gain(3,:) = pow2db((mean(q12(3,:,:),3)+(1-mean(q12(3,:,:),3)).*mean(q11(3,:,:),3)) ./...
                  (mean(q12(1,:,:),3)+(1-mean(q12(1,:,:),3)).*mean(q11(1,:,:),3)));
              
BLER2_gain(1,:) = pow2db(mean(q2(2,:,:),3) ./ mean(q2(1,:,:),3));
BLER2_gain(2,:) = pow2db(mean(q2(3,:,:),3) ./ mean(q2(2,:,:),3));
BLER2_gain(3,:) = pow2db(mean(q2(3,:,:),3) ./ mean(q2(1,:,:),3));


%% Plot

figure (1)

plot(d2-d1,BLER2_gain(1,:),'r', 'linewidth', 1.5);
hold on;grid on;
plot(d2-d1,BLER1_gain(1,:),'b', 'linewidth', 1.5);

plot(d2-d1,BLER2_gain(2,:),'--r', 'linewidth', 1.5);
plot(d2-d1,BLER2_gain(3,:),'-.r', 'linewidth', 1.5);


plot(d2-d1,BLER1_gain(2,:),'--b', 'linewidth', 1.5);
plot(d2-d1,BLER1_gain(3,:),'-.b', 'linewidth', 1.5);

title('BLER gain vs Difference of Distance');
xlabel('Difference of Distance');
ylabel('BLER gain');

legend('Far user','Near user');

figure (2)
semilogy(d2-d1,BLER2_gain(1,:),'r', 'linewidth', 1.5);
hold on;grid on;
semilogy(d2-d1,BLER1_gain(1,:),'b', 'linewidth', 1.5);

semilogy(d2-d1,BLER2_gain(2,:),'--r', 'linewidth', 1.5);
semilogy(d2-d1,BLER2_gain(3,:),'-.r', 'linewidth', 1.5);


semilogy(d2-d1,BLER1_gain(2,:),'--b', 'linewidth', 1.5);
semilogy(d2-d1,BLER1_gain(3,:),'-.b', 'linewidth', 1.5);




toc