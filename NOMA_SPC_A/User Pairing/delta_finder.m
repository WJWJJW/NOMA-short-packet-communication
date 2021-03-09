% This script finds the optimal delta by solving critical point of
% blocklength based on delta

function [opt_delta]=delta_finder(A)
   com_term = (sqrt(3)*sqrt(16*A^3+27*A^2)+9*A)^(1/3);
   term11 = 2^(1/3) * com_term / ((3)^(2/3)*A);
   term12 = 2*2^(2/3) / (3^(1/3) * com_term);
   
   term21 = 2^(1/3) * com_term / (3^(2/3)*A);
   term22 = 4 / (A*sqrt(term11 - term12));
   term23 = term12;
   
   opt_delta = 1/2 * sqrt(-term21 + term22 + term23) - 1/2 * sqrt(term11 - term12);
end