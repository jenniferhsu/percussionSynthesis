function [y, z] = simple_filter(B, A, x, z)
%Y = SIMPLE_FILTER(B, A, x, z)
% this only allows us to do one-pole and one-zero filtering on a vector
% of values (each element of the vector does not have anything to do with
% any other element in the vector)
% this will hopefully make my code run faster...

if length(B) == 1 && length(A) == 2
    % one pole
    y = (B(1)*x - A(2)*z)/A(1);
    z = y;
elseif length(B) == 2 && length(A) == 1
    % one zero
    y = (B(1)*x + B(2)*z)/A(1);
    z = x;
else
    disp('I can only do one-pole and one-zero filters...sorry');
end


end