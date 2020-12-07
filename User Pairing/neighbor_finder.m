function [neighbor_1, neighbor_2] = neighbor_finder(cur_combinition, K)
    
    while 1
        pair_idx = randi(K,1,2);
        if pair_idx(1) ~= pair_idx(2)
            break;
        end
    end
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
    % output
    neighbor_1 = cur_combinition;
    neighbor_1(pair_idx(1),:) = pairA;
    neighbor_1(pair_idx(2),:) = pairB;
    
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
    % output
    neighbor_2 = cur_combinition;
    neighbor_2(pair_idx(1),:) = pairA;
    neighbor_2(pair_idx(2),:) = pairB;
end