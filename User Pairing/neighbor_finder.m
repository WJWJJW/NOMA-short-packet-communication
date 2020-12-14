function [neighbor_1, neighbor_2, pair_idx] = neighbor_finder(cur_combinition, K)
    
    % random pick
    pair_idx = zeros(1,2,'uint32');
    pair_idx(1) = randi(K,1);
    valid_choices = setdiff (1:K, pair_idx(1));
    pair_idx(2) = valid_choices(randi(K-1));

    
    % neighbor 1
    pairA = cur_combinition(pair_idx(1),:);
    pairB = cur_combinition(pair_idx(2),:);
    % swap
    tmp_user = pairA(1,1);
    pairA(1,1) = pairB(1,1);
    pairB(1,1) = tmp_user;
    % sort again
    pairA = sort(pairA);
    pairB = sort(pairB);
%     % output
%     neighbor_1 = cur_combinition;
%     neighbor_1(pair_idx(1),:) = pairA;
%     neighbor_1(pair_idx(2),:) = pairB;
    neighbor_1 = [pairA; pairB];
    
    % neighbor 2
    pairA = cur_combinition(pair_idx(1),:);
    pairB = cur_combinition(pair_idx(2),:);
    % swap
    tmp_user = pairA(1,2);
    pairA(1,2) = pairB(1,1);
    pairB(1,1) = tmp_user;
    % sort again
    pairA = sort(pairA);
    pairB = sort(pairB);
%     % output
%     neighbor_2 = cur_combinition;
%     neighbor_2(pair_idx(1),:) = pairA;
%     neighbor_2(pair_idx(2),:) = pairB;
    neighbor_2 = [pairA; pairB];
end