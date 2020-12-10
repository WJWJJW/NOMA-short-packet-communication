clc; clear variables; close all;
N = 1e6; % number of Monte Carlo
K = 5;  % number of cluster (number of user  = 2K)
NN = 80; % number of information bit
N1 = NN;
N2 = NN;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = -10:2:0;                    %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt./ no;
RHO = pow2db(rho);

eta = 4;
dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);

user_distance = randi([10 330],1,2*K);
user_distance = sort(user_distance);

RP_user_pairing = zeros(K,2);
User_pre_grouping = zeros(K,2);
Timer = zeros(2, length(Pt));

for u=1:length(Pt)
    %% Random Paring (RP)
    RP_indices = randperm(2*K);
    for ii=1:K
        RP_user_pairing(ii,:,u) = sort(user_distance(RP_indices(2*ii-1:2*ii)));
    end
    
    [RP_thred1_check(:,u), RP_thred2_check(:,u)] = thred_checker(RP_user_pairing(:,:,u), K, eplsion1R, eplsion2R, rho(u), eta);
    
    % Total blocklength for random pairing
    [sum_RP_opt_M(u), RP_opt_M(:,u)] = M_cal(N1, RP_user_pairing(:,:,u), K,...
                                        eplsion1R,eplsion2R,rho(u),N,eta);
    
    %% User Pre-Grouping
    for ii=1:K
        User_pre_grouping(ii,1,u) = user_distance(ii);
        User_pre_grouping(ii,2,u) = user_distance(K -1 + ii);
    end
    
    [UPG_thred1_check(:,u), UPG_thred2_check(:,u)] = thred_checker(User_pre_grouping(:,:,u), K, eplsion1R, eplsion2R, rho(u), eta);
    
    % Total blocklength for User Pre-Grouping
    [sum_UPG_opt_M(u), UPG_opt_M(:,u)] = M_cal(N1, User_pre_grouping(:,:,u), K,...
                                        eplsion1R,eplsion2R,rho(u),N,eta);
    
    tic;
    %% Hill Climbing Pairing
    % RP and calculate blocklength as current optimum blocklength
    [sum_HCP_opt_M(u), HCP_opt_M(:,u)] = M_cal(N1, RP_user_pairing(:,:,u), K,...
                                        eplsion1R,eplsion2R,rho(u),N,eta);
    
    cur_combinition = RP_user_pairing(:,:,u); 
    while 1
        % find neighbor
        [neighbor_1, neighbor_2] = neighbor_finder(cur_combinition, K);
        [sum_nei1_opt_M, nei1_opt_M(:,u)] = M_cal(N1, neighbor_1, K,...
                                        eplsion1R,eplsion2R,rho(u),N,eta);
        [sum_nei2_opt_M, nei2_opt_M(:,u)] = M_cal(N1, neighbor_2, K,...
                                        eplsion1R,eplsion2R,rho(u),N,eta);
        
        % Find the best neighbor
        if sum_nei1_opt_M < sum_nei2_opt_M
            sum_nei_opt_M = sum_nei1_opt_M;
            nei_opt_M = nei1_opt_M(:,u);
            neighbor = neighbor_1;
        else
            sum_nei_opt_M = sum_nei2_opt_M;
            nei_opt_M = nei2_opt_M(:,u);
            neighbor = neighbor_2;
        end
        
        % Find the solution for this iteration

        if  sum_nei_opt_M <= sum_HCP_opt_M (u)
            sum_HCP_opt_M (u) = sum_nei_opt_M;
            HCP_opt_M(:,u) = nei_opt_M;
            cur_combinition = neighbor;
        else
            break;
        end
        
            
    end % End HCP
    toc;
    Timer(1,u) = toc;
    
    tic;
    %% Simulated Annealing Pairing
    % Initialization
    Temperature = 20;
    Temperature_min = 0.1;
    annealing_factor = 0.7; 
    Time_budget = 50;
    cur_time = 0;
    % RP and calculate blocklength as current optimum blocklength
    [sum_SAP_opt_M(u), SAP_opt_M(:,u)] = M_cal(N1, RP_user_pairing(:,:,u), K,...
                                        eplsion1R,eplsion2R,rho(u),N,eta);
                                    
    cur_combinition = RP_user_pairing(:,:,u); 
    while 1
        % Time update
        cur_time = cur_time+1;
        % Time budget check
        if cur_time > Time_budget
            break;
        end
        % Find neighbor
        [neighbor_1, neighbor_2] = neighbor_finder(cur_combinition, K);
        % Choose neighbor randomly
        if randi(2) == 1
            neighbor = neighbor_1;
        else
            neighbor = neighbor_2;
        end
        % Calculate blocklength for neighbor
        [sum_nei_opt_M, nei_opt_M] = M_cal(N1, neighbor, K,...
                                        eplsion1R,eplsion2R,rho(u),N,eta);
        
        % Find the solution for this iteration
        % Better solution
        if  sum_nei_opt_M <= sum_SAP_opt_M (u)
            sum_SAP_opt_M (u) = sum_nei_opt_M;
            SAP_opt_M(:,u) = nei_opt_M;
            cur_combinition = neighbor;
        % Worse solution
        else
            delta_E = sum_SAP_opt_M (u) - sum_nei_opt_M;
            rn = rand(1);
            % Accept worse solution
            if exp(delta_E/Temperature) >= rn
                sum_SAP_opt_M (u) = sum_nei_opt_M;
                SAP_opt_M(:,u) = nei_opt_M;
                cur_combinition = neighbor;
            end           
        end
        % Annealing
        Temperature = annealing_factor * Temperature;
        if Temperature < Temperature_min
            break;
        end
            
    end % End SAP
    Simulated_Anealing_Pairing(:,:,u) = cur_combinition;
     [SAG_thred1_check(:,u), SAG_thred2_check(:,u)] = thred_checker(Simulated_Anealing_Pairing(:,:,u), K, eplsion1R, eplsion2R, rho(u), eta);
    toc;
    Timer(2,u) = toc;
end

figure (1)

plot(Pt, sum_RP_opt_M,'b');
hold on; grid on;
plot(Pt, sum_UPG_opt_M,'g');
plot(Pt, sum_HCP_opt_M,'r');
plot(Pt, sum_SAP_opt_M,'Color',[1 0.5 0]);


ylabel('blocklength');
legend('Random Pairing','User Pre-Grouping',...
       'Hill Climbing Based Pairing', 'Simulated Anealing Based Pairing');