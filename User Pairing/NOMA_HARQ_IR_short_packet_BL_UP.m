clc; clear variables; close all;

N = 1e6; % number of channel tap
NNN = 2; % number of Monte Carlo
K = 5;  % number of cluster (number of user  = 2K)
NN = 80; % number of information bit
N1 = NN;
N2 = NN;

eplsion1R = 10^-5;
eplsion2R = 10^-4;

Pt = -10:2:10;                    %Transmit Power in dBm
pt = (10^-3).*db2pow(Pt);    %Transmit Power (linear scale)

BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)

rho = pt./ no;
RHO = pow2db(rho);

eta = 4;
dis_thred1 = (eplsion2R/eplsion1R)^(1/eta);

RP_user_pairing = zeros(K,2,length(Pt));
User_pre_grouping = zeros(K,2,length(Pt));
Simulated_Anealing_Pairing = zeros(K,2,length(Pt));

sum_RP_opt_M_j = zeros(NNN,length(Pt));
sum_UPG_opt_M_j = zeros(NNN,length(Pt));
sum_HCP_opt_M_j = zeros(NNN,length(Pt));
sum_SAP_opt_M_j = zeros(NNN,length(Pt));

RP_opt_M = zeros(K);
HCP_opt_M = zeros(K);
SAP_opt_M = zeros(K);

Timer_RP = zeros(NNN,length(Pt));
Timer_UPG = zeros(NNN,length(Pt));
Timer_HCP = zeros(NNN,length(Pt));
Timer_SAP = zeros(NNN,length(Pt));


s = tic;
parfor u=1:length(Pt)
    h = (randn(1,N)+1i*randn(1,N));
    lamda = mean(abs(h).^2);
    for jj = 1:NNN
        % Generate user randomly
        user_distance = randi([10 330],1,2*K);
        user_distance = sort(user_distance);
        
        a = tic;
        %% Random Paring (RP)
        RP_indices = randperm(2*K);
        tmp = zeros(K,2);
        for ii=1:K
            tmp(ii,:) = sort(user_distance(RP_indices(2*ii-1:2*ii)));
%             RP_user_pairing_thred_check(ii,u) = double((RP_user_pairing(ii,2,u) / RP_user_pairing(ii,1,u)) > dis_thred1);
        end
        RP_user_pairing(:,:,u) = tmp;
        % Total blocklength for random pairing
        [sum_RP_opt_M_j(jj,u), RP_opt_M] = M_cal(N1, RP_user_pairing(:,:,u), K,...
                                        eplsion1R,eplsion2R,rho(u),eta,lamda); 

        Timer_RP(jj,u) = toc(a);                               
                         
        b = tic;       
        %% User Pre-Grouping
        tmp = zeros(K,2);
        for ii=1:K
            tmp(ii,1) = user_distance(ii);
            tmp(ii,2) = user_distance(K -1 + ii);
%             User_pre_grouping_thred_check(ii,u) = double((User_pre_grouping(ii,2,u) / User_pre_grouping(ii,1,u)) > dis_thred1);
        end
        
        User_pre_grouping(:,:,u) = tmp;
        
        % Total blocklength for User Pre-Grouping
        [sum_UPG_opt_M_j(jj,u)] = M_cal(N1, User_pre_grouping(:,:,u), K,...
                                            eplsion1R,eplsion2R,rho(u),eta,lamda);                         
        
        Timer_UPG(jj,u) = toc(b);
        
        c = tic;
        %% Hill Climbing Pairing
        % RP and calculate blocklength as current optimum blocklength
        sum_HCP_opt_M_j(jj,u) = sum_RP_opt_M_j(jj,u);
        HCP_opt_M = RP_opt_M;

        cur_combinition = RP_user_pairing(:,:,u); 
        while 1
            % find neighbor
            [neighbor_1, neighbor_2, diff_idx] = neighbor_finder(cur_combinition, K);
            % calculate sum of non-changing pair
            tmp_sum = sum_HCP_opt_M_j(jj,u) - HCP_opt_M(diff_idx(1)) - HCP_opt_M(diff_idx(2));

            % calculate sum of changing pair
            [sum_nei1_opt_M, nei1_opt_M] = M_cal(N1, neighbor_1, 2,...
                                            eplsion1R,eplsion2R,rho(u),eta,lamda);
            sum_nei1_opt_M = tmp_sum + sum_nei1_opt_M;


            [sum_nei2_opt_M, nei2_opt_M] = M_cal(N1, neighbor_2, 2,...
                                            eplsion1R,eplsion2R,rho(u),eta,lamda);
            sum_nei2_opt_M = tmp_sum + sum_nei2_opt_M;


            % Find the best neighbor
            if sum_nei1_opt_M < sum_nei2_opt_M
                sum_nei_opt_M = sum_nei1_opt_M;
                nei_opt_M = nei1_opt_M;
                neighbor = neighbor_1;
            else
                sum_nei_opt_M = sum_nei2_opt_M;
                nei_opt_M = nei2_opt_M;
                neighbor = neighbor_2;
            end

            % Find the solution for this iteration

            if  sum_nei_opt_M < sum_HCP_opt_M_j (jj,u)
                sum_HCP_opt_M_j (jj,u) = sum_nei_opt_M;
                HCP_opt_M(diff_idx(1)) = nei_opt_M(1);
                HCP_opt_M(diff_idx(2)) = nei_opt_M(2);
                cur_combinition(diff_idx(1),:) = neighbor(1,:);
                cur_combinition(diff_idx(2),:) = neighbor(2,:);
            else
                break;
            end
        end % End HCP

        
        Timer_HCP(jj,u) = toc(c);
    
        d = tic;
        %% Simulated Annealing Pairing
        % Initialization
        Temperature = 20;
        Temperature_min = 0.1;
        annealing_factor = 0.7; 
        Time_budget = 50;
        cur_time = 0;
        % RP and calculate blocklength as current optimum blocklength
        sum_SAP_opt_M_j(jj,u) = sum_RP_opt_M_j(jj,u);
        SAP_opt_M = RP_opt_M;

        cur_combinition = RP_user_pairing(:,:,u); 

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
            tmp_sum = sum_SAP_opt_M_j(jj,u) - SAP_opt_M(diff_idx(1)) - SAP_opt_M(diff_idx(2));

            % calculate sum of changing pair
            [sum_nei_opt_M, nei_opt_M] = M_cal(N1, neighbor, 2,...
                                            eplsion1R,eplsion2R,rho(u),eta,lamda);
            sum_nei_opt_M = tmp_sum + sum_nei_opt_M;

            % Find the solution for this iteration
            % Better solution
            if  sum_nei_opt_M <= sum_SAP_opt_M_j (jj,u)
                sum_SAP_opt_M_j (jj,u) = sum_nei_opt_M;
                SAP_opt_M(diff_idx(1)) = nei_opt_M(1);
                SAP_opt_M(diff_idx(2)) = nei_opt_M(2);
                cur_combinition(diff_idx(1),:) = neighbor(1,:);
                cur_combinition(diff_idx(2),:) = neighbor(2,:);
            % Worse solution
            else
                delta_E = sum_SAP_opt_M_j (jj,u) - sum_nei_opt_M;
                rn = rand(1);
                % Accept worse solution
                if exp(delta_E/Temperature) >= rn
                    sum_SAP_opt_M_j (jj,u) = sum_nei_opt_M;
                    SAP_opt_M(diff_idx(1)) = nei_opt_M(1);
                    SAP_opt_M(diff_idx(2)) = nei_opt_M(2);
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
        
        Timer_SAP(jj,u) = toc(d);
    Simulated_Anealing_Pairing(:,:,u) = cur_combinition;                                                                        
    end

end
sum_RP_opt_M = mean(sum_RP_opt_M_j);
sum_UPG_opt_M = mean(sum_UPG_opt_M_j);
sum_HCP_opt_M = mean(sum_HCP_opt_M_j);
sum_SAP_opt_M = mean(sum_SAP_opt_M_j);

t = toc(s)


figure (1)

plot(Pt, sum_RP_opt_M,'b');
hold on; grid on;
plot(Pt, sum_UPG_opt_M,'g');
plot(Pt, sum_HCP_opt_M,'r');
plot(Pt, sum_SAP_opt_M,'Color',[1 0.5 0]);


ylabel('blocklength');
legend('Random Pairing','User Pre-Grouping','Hill Climbing Based Pairing','Simulated Anealing Based Pairing');

name_str = ['UP_test_' num2str(NNN) 'times.png'];
saveas(gcf,name_str);

saveas(gcf,'UP_test.png');