function f = membraneRectModes(dx, dy, T, sig)
%%f = membraneRectModes(dx, dy, T, sig)
% calculates the modal frequencies of a rectangular membrane of size (dx, dy)
% with tension T and mass density per unit area sig. The rectangular membrane
% is assumed to have inverting reflection on all four sides.
% The equation for calculating the modal frequencies is from Morse
% and Ingard's Theoretical Acoustics (1968) on page 205.
%
% example:
% f = membraneRectModes (3, 3, 1, 1);

% total number of modes in each direction (number of half waves)
% this is just what makes sense for now. i might change this depending on
% resulting sounds
Nm = floor(dx/2) + 1;
Nn = floor(dy/2) + 1;

% store the modal frequencies in f
f = zeros(Nm, Nn);

for m=1:Nm
    for n=1:Nn
        f(m,n) = 0.5 * sqrt(T/sig) * sqrt((m^2/dx^2) + (n^2/dy^2));
    end
end

% normalize frequencies
f = f/f(1);

end


