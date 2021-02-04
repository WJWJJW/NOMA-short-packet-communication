function [i, j] = linear_sum_assignment(cost_matrix)
% function [i, j] = linear_sum_assignment(cost_matrix)
%
% Hungarian algorithm (Kuhn-Munkres) for solving the linear sum assignment
% problem.
%
% Input: 
% COST_MATRIX: n-by-m matrix with costs. Rectangular cost matrices 
%              are supported. 
% 
% Output: 
% I, J: matching.
% Row I(k) is matched to column J(k), 
%     k = 1, 2, 3, ... min(size(COST_MATRIX)). 
%

% Based on python implementation by Brian Clapper and Gael Varoquaux.
% Adapted to matlab by Ondrej Drbohlav.
%
% Copyright (c) 2008 Brian M. Clapper <bmc@clapper.org>, Gael Varoquaux
% Copyright (c) 2016 Ondrej Drbohlav <drbohlav@fel.cvut.cz>
% License: 3-clause BSD

True = 1;
False = 0;

if length(size(cost_matrix)) ~= 2
    error(sprintf('Expected a matrix (2-d array), got a %i-d array'), ...
          length(size(cost_matrix))); 
end

%  The algorithm expects more columns than rows in the cost matrix.
if size(cost_matrix, 2) < size(cost_matrix, 1)
    cost_matrix = cost_matrix';
    transposed = True;
else
    transposed = False;
end

state = Hungary_(cost_matrix); 

step = @step1; 

while ~isequal(step, @none)
    [step, state] = step(state); 
end

if transposed
    marked = state.marked';
else
    marked = state.marked;
end
    
[i, j] = find(marked == 1); 


function s = Hungary_(cost_matrix)
[n, m] = size(cost_matrix); 

s = struct('C', cost_matrix, ...
           'row_uncovered', true(n, 1), ...
           'col_uncovered', true(1, m), ...
           'Z0_r', 1, ... 
           'Z0_c', 1, ... 
           'path', ones(n+m, 2, 'int64'), ...
           'marked', zeros(n, m, 'int8')); 

function none()
return 

function state = clear_covers(state)
True = 1;
False = 0;
% """Clear all covered matrix cells"""
state.row_uncovered(:) = True; 
state.col_uncovered(:) = True;

%# Individual steps of the algorithm follow, as a state machine: they return
%# the next step to be taken (function to be called), if any.

function [fnhandle, state] = step1(state)
True = 1;
False = 0; 
%"""Steps 1 and 2 in the Wikipedia page."""
%
%    # Step 1: For each row of the matrix, find the smallest element and
%    # subtract it from every element in its row.
state.C = state.C - ...
          repmat(min(state.C,[],2), [1, size(state.C, 2)]); 
%    # Step 2: Find a zero (Z) in the resulting matrix. If there is no
%    # starred zero in its row or column, star Z. Repeat for each element
%    # in the matrix.
[i_, j_] = find(state.C == 0);

for k = 1:length(i_)
    i = i_(k); 
    j = j_(k); 
    if state.col_uncovered(j) & state.row_uncovered(i)
        state.marked(i, j) = 1;
        state.col_uncovered(j) = False;
        state.row_uncovered(i) = False;
    end
end

state = clear_covers(state);
fnhandle = @step3;

function [fnhandle, state] = step3(state)
True = 1;
False = 0;
%    """
%    Cover each column containing a starred zero. If n columns are covered,
%    the starred zeros describe a complete set of unique assignments.
%    In this case, Go to DONE, otherwise, Go to Step 4.
%    """
marked = (state.marked == 1);
state.col_uncovered( any(marked) ) = False; 

if sum(marked(:)) < size(state.C, 1)
    fnhandle = @step4; 
else
    fnhandle = @none; 
end
        
function [fnhandle, state] = step4(state)
True = 1;
False = 0;
%    """
%    Find a noncovered zero and prime it. If there is no starred zero
%    in the row containing this primed zero, Go to Step 5. Otherwise,
%    cover this row and uncover the column containing the starred
%    zero. Continue in this manner until there are no uncovered zeros
%    left. Save the smallest uncovered value and Go to Step 6.
%    """
%    # We convert to int as numpy operations are faster on int

[n, m] = size(state.C); 

% covered_C: bool
covered_C = logical(state.C == 0);
covered_C(~state.row_uncovered, :) = 0;
covered_C(:, ~state.col_uncovered, :) = 0;

while True
    %# Find an uncovered zero
    idx_zero = find(covered_C, 1); 
    
    if length(idx_zero) == 0
        fnhandle = @step6; 
        return 
    else
        row = mod(idx_zero-1, n)+1; col = floor((idx_zero-1)/n)+1; 
        state.marked(row, col) = 2;
        %# Find the first starred element in the row
        star_col = find(state.marked(row, :) == 1, 1); 
        if isempty(star_col)
            %# Could not find one
            state.Z0_r = row;
            state.Z0_c = col;
            fnhandle = @step5; 
            return
        else
            col = star_col;
            state.row_uncovered(row) = False;
            state.col_uncovered(col) = True;
            covered_C(~state.row_uncovered, col) = 0 ; 
            covered_C(row, :) = 0; 
        end
    end
end


function [fnhandle, state] = step5(state)
True = 1;
False = 0;
%    """
%    Construct a series of alternating primed and starred zeros as follows.
%    Let Z0 represent the uncovered primed zero found in Step 4.
%    Let Z1 denote the starred zero in the column of Z0 (if any).
%    Let Z2 denote the primed zero in the row of Z1 (there will always be one%).
%    Continue until the series terminates at a primed zero that has no starred
%    zero in its column. Unstar each starred zero of the series, star each
%    primed zero of the series, erase all primes and uncover every line in th%e
%    matrix. Return to Step 3
%    """
el0 = 1;
el1 = 2;
[n, m] = size(state.C); 

count = 1; 
state
state.path(count, el0) = state.Z0_r;
state.path(count, el1) = state.Z0_c;

while True
    %# Find the first starred element in the col defined by
    %# the path.
    row = find(state.marked(:, state.path(count, el1) ) == 1, 1); 
    if isempty(row)
        %# Could not find one
        break; 
    else
        count = count + 1; 
        state.path(count, el0) = row; 
        state.path(count, el1) = state.path(count - 1, el1);
    end
        
    %# Find the first prime element in the row defined by the
    %# first path step
    col = find(state.marked( state.path(count, el0), : ) == 2, 1); 
    % Note: there will always be such col, thus we do not test 
    % for isempty(col)
    count = count + 1; 
    state.path(count, el0) = state.path(count - 1, el0);
    state.path(count, el1) = col;
end
    
%# Convert paths:
for i = 1:(count)
    if state.marked( state.path(i, el0), state.path(i, el1) ) == 1
        state.marked( state.path(i, el0), state.path(i, el1)) = 0; 
    else
        state.marked( state.path(i, el0), state.path(i, el1) ) = 1;
    end
end
    
state = clear_covers(state);
%# Erase all prime markings
state.marked(state.marked == 2) = 0;
fnhandle = @step3; return 


function [fnhandle, state] = step6(state)
%    """
%    Add the value found in Step 4 to every element of each covered row,
%    and subtract it from every element of each uncovered column.
%    Return to Step 4 without altering any stars, primes, or covered lines.
%    """
%    # the smallest uncovered value in the matrix
if any(state.row_uncovered) & any(state.col_uncovered)
    minval = min(min(state.C(state.row_uncovered, ...
                             state.col_uncovered))); 
    state.C(~state.row_uncovered, :) =  state.C(~state.row_uncovered, :) ...
        + minval;
    state.C(:, state.col_uncovered) = state.C(:, state.col_uncovered) ...
        - minval; 
end

fnhandle = @step4;