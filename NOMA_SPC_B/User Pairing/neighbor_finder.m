function [neighbor_1_idx, neighbor_2_idx] = neighbor_finder(cur_combinition_idx, K)
    
    % random pick
    pair_idx = zeros(1,2,'uint32');
    pair_idx(1) = randi(K,1);
    valid_choices = setdiff (1:K, pair_idx(1));
    pair_idx(2) = valid_choices(randi(K-1));

    % swap by linear indexing method
    % Neighbor 1
    neighbor_1_idx = cur_combinition_idx;
    neighbor_1_idx([pair_idx(1) pair_idx(2)]) ...
        = neighbor_1_idx([pair_idx(2) pair_idx(1)]);
    % Neighbor 2
    neighbor_2_idx = cur_combinition_idx;
    neighbor_2_idx([pair_idx(1) pair_idx(2)+K]) ...
        = neighbor_2_idx([pair_idx(2)+K pair_idx(1)]);
end