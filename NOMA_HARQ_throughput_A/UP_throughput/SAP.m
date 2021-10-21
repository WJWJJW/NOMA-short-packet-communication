% Simulated Annealing Pairing (SAP)
function [sum_opt_M, opt_M, Simulated_Anealing_Pairing, cur_combinition_idx]=SAP(user_distance, N, K, target_BLER, rho, eta, lamda)
    % Initialization
    Temperature = 50;
    Temperature_min = 0.0001;
    annealing_factor = 0.9; 
    Time_budget = 100000;
    cur_time = 0;
    % RP and calculate blocklength as current optimum blocklength
    RP_indices = randperm(2*K);
    RP_indices = sort(reshape(RP_indices,K,2),2);
    RP_user_pairing = user_distance(RP_indices);
    target_BLER_pair = target_BLER(RP_indices);

    
    % Total blocklength for random pairing
    [sum_opt_M, opt_M] = M_cal_Mod(N,RP_user_pairing, K, target_BLER_pair,rho,eta,lamda);
    
    cur_combinition_idx = RP_indices;
                                    
    while 1
        % Time update
        cur_time = cur_time+1;
        % Time budget check
        if cur_time > Time_budget
            break;
        end
        % Find neighbor
        [neighbor_1_idx, neighbor_2_idx] = neighbor_finder(cur_combinition_idx, K);
        % Choose neighbor randomly
        if randi(2) == 1
            neighbor_idx = neighbor_1_idx;
        else
            neighbor_idx = neighbor_2_idx;
        end

        % calculate blocklength of new combination
        [sum_nei_opt_M, nei_opt_M] = ...
               M_cal_Mod(N, user_distance(neighbor_idx), K,target_BLER(neighbor_idx),rho,eta,lamda);
        
        % Find the solution for this iteration
        % Better solution
        if  sum_nei_opt_M < sum_opt_M
            sum_opt_M = sum_nei_opt_M;
            opt_M = nei_opt_M;
            cur_combinition_idx = neighbor_idx;

        % Worse solution
        else
            delta_E = sum_opt_M - sum_nei_opt_M;
            rn = rand(1);
            % Accept worse solution
            if exp(delta_E/Temperature) >= rn
                sum_opt_M = sum_nei_opt_M;
                opt_M = nei_opt_M;
                cur_combinition_idx = neighbor_idx;
            end           
        end
        % Annealing
        Temperature = annealing_factor * Temperature;
        if Temperature < Temperature_min
            break;
        end
            
    end % End SAP
    Simulated_Anealing_Pairing = user_distance(cur_combinition_idx);
end