
% Transmitted Power
Pt = 20:10:30;               %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

% AWGN
% BW = 10^7;                  %System bandwidth
% No = -174 + 10*log10(BW);   %Noise power (dBm)
% no = (10^-3)*10.^(No/10);   %Noise power (linear scale)
No = -100;
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt./ no;
RHO = pow2db(rho);

beta = 0.5;
OMA_PA = 0.5;

eta = 4;
target_BLER = [1e-5 1e-4];

user_distance = 50:2:250;

d_min = zeros(4,length(user_distance),length(Pt));

for u = 1:length(user_distance)
    d_min(1,u,1) = min_d_finder(user_distance(u), target_BLER(1), target_BLER(2), rho(1));
    d_min(2,u,1) = min_d_finder(user_distance(u), target_BLER(1), target_BLER(1), rho(1));
    d_min(3,u,1) = min_d_finder(user_distance(u), target_BLER(2), target_BLER(2), rho(1));
    d_min(4,u,1) = min_d_finder(user_distance(u), target_BLER(2), target_BLER(1), rho(1));

    d_min(1,u,2) = min_d_finder(user_distance(u), target_BLER(1), target_BLER(2), rho(2));
    d_min(2,u,2) = min_d_finder(user_distance(u), target_BLER(1), target_BLER(1), rho(2));
    d_min(3,u,2) = min_d_finder(user_distance(u), target_BLER(2), target_BLER(2), rho(2));
    d_min(4,u,2) = min_d_finder(user_distance(u), target_BLER(2), target_BLER(1), rho(2));

end

figure(1)
plot(user_distance, d_min(1,:,1),'b');
hold on;grid on;
plot(user_distance, d_min(2,:,1),'r');
plot(user_distance, d_min(3,:,1),'g');
plot(user_distance, d_min(4,:,1),'c');

plot(user_distance, d_min(1,:,2),'--b');
plot(user_distance, d_min(2,:,2),'--r');
plot(user_distance, d_min(3,:,2),'--g');
plot(user_distance, d_min(4,:,2),'--c');

xlabel('NU distance');
ylabel('Min. distance for FU');
legend('\epsilon_i = 10^{-5}, \epsilon_j = 10^{-4}',...
       '\epsilon_i = 10^{-5}, \epsilon_j = 10^{-5}',...
       '\epsilon_i = 10^{-4}, \epsilon_j = 10^{-4}',...
       '\epsilon_i = 10^{-4}, \epsilon_j = 10^{-5}');
set(gca, 'FontName', 'Times New Roman'); 

