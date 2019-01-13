function [B,A] = find_onepoleBA_given_z(p)
% given a pole location p for a one pole filter, find the B and A filter
% coefficients that will keep the final gain at 1 or below

%     a1 = -1/z;
%     w = 0:0.01:2*pi;    % vector of frequencies
%     b0 = sqrt(1 + a1^2 + 2*a1*(cos(w)));
%     
%     % plot(w, b0);
%     
%     B = -min(b0);
%     A = [1 a1];
%     %figure;
%     %freqz(B, A, 1024)


% let's try the Stk version
    if p > 0
        b0 = 1-p;
    else
        b0 = 1+p;
    end
    a1 = -p;
    
    B = b0;
    A = [1, a1];


end