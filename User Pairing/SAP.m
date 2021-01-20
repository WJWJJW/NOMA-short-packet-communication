% Simulated Annealing Pairing (SAP)
function [sum_opt_M, opt_M, Simulated_Anealing_Pairing]=SAP(user_distance, N, K, eplsion1R, eplsion2R, rho, eta, lamda, delta)
    % Initialization
    Temperature = 50;
    Temperature_min = 0.0001;
    annealing_factor = 0.9; 
    Time_budget = 100;
    cur_time = 0;
    % RP and calculate blocklength as current optimum blocklength
    [sum_opt_M, opt_M, cur_combinition] = RP(user_distance, N, K, eplsion1R, eplsion2R, rho, eta, lamda, delta);
                                    
    while 1
        % Time update
        cur_time = cur_time+1;
        % Time budget check
        if cur_time > Time_budget
            break;
        end
        % Find neighbor
        [neighbor_1, neighbor_2, diff_idx] = neighbor_finder(cur_combinition, K);
        % Choose neighbor randomly
        if randi(2) == 1
            neighbor = neighbor_1;
        else
            neighbor = neighbor_2;
        end

        % calculate sum of non-changing pair
        tmp_sum = sum_opt_M - opt_M(diff_idx(1)) - opt_M(diff_idx(2));
        
        % calculate sum of changing pair
        [sum_nei_opt_M, nei_opt_M] = M_cal(N, neighbor, 2,eplsion1R,eplsion2R,rho,eta,lamda, delta);
        sum_nei_opt_M = tmp_sum + sum_nei_opt_M;
        
        % Find the solution for this iteration
        % Better solution
        if  sum_nei_opt_M < sum_opt_M
            sum_opt_M = sum_nei_opt_M;
            opt_M(diff_idx(1)) = nei_opt_M(1);
            opt_M(diff_idx(2)) = nei_opt_M(2);
            cur_combinition(diff_idx(1),:) = neighbor(1,:);
            cur_combinition(diff_idx(2),:) = neighbor(2,:);
        % Worse solution
        else
            delta_E = sum_opt_M - sum_nei_opt_M;
            rn = rand(1);
            % Accept worse solution
            if exp(delta_E/Temperature) >= rn
                sum_opt_M = sum_nei_opt_M;
                opt_M(diff_idx(1)) = nei_opt_M(1);
                opt_M(diff_idx(2)) = nei_opt_M(2);
                cur_combinition(diff_idx(1),:) = neighbor(1,:);
                cur_combinition(diff_idx(2),:) = neighbor(2,:);
            end           
        end
        % Annealing
        Temperature = annealing_factor * Temperature;
        if Temperature < Temperature_min
            break;
        end
            
    end % End SAP
    Simulated_Anealing_Pairing = cur_combinition;
end