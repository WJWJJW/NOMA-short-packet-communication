clc; clear variables; close all;
N = 1e6; % number of Monte Carlo
ncluster = 5:10;  % number of cluster (number of user  = 2K)
NN = 80; % number of information bit
N1 = NN;
N2 = NN;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = 0;                    %Transmit Power in dBm
pt = (10^-3)*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10^(No/10);   %Noise power (linear scale)

rho = pt/ no;
RHO = pow2db(rho);

eta = 4;

for u=1:length(ncluster)
    K = ncluster(u);
    if u == 1
        user_distance = randi([10 330],1,2*K);
    else
        user_distance_add = randi([10 330],1,2);
        user_distance = [user_distance user_distance_add];
    end
    
    user_distance = sort(user_distance)
    RP_user_pairing = zeros(K,2);
    User_pre_grouping = zeros(K,2);
    
    %% Random Paring (RP)
    RP_indices = randperm(2*K);
    for ii=1:K
        RP_user_pairing(ii,:) = sort(user_distance(RP_indices(2*ii-1:2*ii)));
    end
    
    %% User Pre-Grouping
    for ii=1:K
        User_pre_grouping(ii,1) = user_distance(ii);
        User_pre_grouping(ii,2) = user_distance(K -1 + ii);
    end
    
    % Total blocklength for User Pre-Grouping
    [sum_UPG_opt_M(u)] = M_cal(N1, User_pre_grouping(:,:), K,...
                                        eplsion1R,eplsion2R,rho,N,eta);
    
    %% Simulated Annealing Pairing
    % Initialization
    Temperature = 20;
    Temperature_min = 0.1;
    annealing_factor = 0.7; 
    Time_budget = 50;
    cur_time = 0;
    % RP and calculate blocklength as current optimum blocklength
    [sum_SAP_opt_M(u)] = M_cal(N1, RP_user_pairing(:,:), K,...
                                        eplsion1R,eplsion2R,rho,N,eta);
                                    
    cur_combinition = RP_user_pairing(:,:); 
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
        [sum_nei_opt_M] = M_cal(N1, neighbor, K,...
                                        eplsion1R,eplsion2R,rho,N,eta);
        
        % Find the solution for this iteration
        % Better solution
        if  sum_nei_opt_M <= sum_SAP_opt_M (u)
            sum_SAP_opt_M (u) = sum_nei_opt_M;
            cur_combinition = neighbor;
        % Worse solution
        else
            delta_E = sum_SAP_opt_M (u) - sum_nei_opt_M;
            rn = rand(1);
            % Accept worse solution
            if exp(delta_E/Temperature) >= rn
                sum_SAP_opt_M (u) = sum_nei_opt_M;
                cur_combinition = neighbor;
            end           
        end
        % Annealing
        Temperature = annealing_factor * Temperature;
        if Temperature < Temperature_min
            break;
        end
            
    end % End SAP
    
    
end

figure (1)

plot(ncluster, sum_UPG_opt_M,'g');
hold on; grid on;
plot(ncluster, sum_SAP_opt_M,'Color',[1 0.5 0]);


ylabel('blocklength');
legend('User Pre-Grouping','Simulated Anealing Based Pairing');