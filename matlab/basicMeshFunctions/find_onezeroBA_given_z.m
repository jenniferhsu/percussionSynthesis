function [B,A] = find_onezeroBA_given_z(z)
% given a pole location z for a one zero filter, find the B and A filter
% coefficients that will keep the final gain at 1 or below
% see the Stk code for oneZero.setZero()

% you can figure this out from the gain formula for the one-zero filter on
% julius's website: 
%   https://ccrma.stanford.edu/~jos/fp/One_Zero.html
% and then finding two equations in terms of b1 and z (make b0 in terms of
% b1 and z) and then substituting
%   cos(wT)  for max(cos(wT)) = 1

    if z > 0
        b0 = 1/(1+z);
    else
        b0 = 1/(1-z);
    end
    b1 = -b0*z;
    
    B = [b0 b1];
    A = 1;
    %figure;
    %freqz(B, A, 1024)

end
