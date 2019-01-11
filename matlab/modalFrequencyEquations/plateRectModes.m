function f = plateRectModes(dx, dy, h, cL)
%%f = plateRectModes (dx, dy, h, cL)
% calculates the modal frequencies of a rectangular plate of size (dx, dy)
% with mesh thickness h and longitudinal wave speed cL. The rectangular plate
% is assumed to be simply-supported (non-inverting reflection on all four sides).
% The equation for calculating the modal frequencies is from Tom
% Rossing's The Science of Percussion Instruments (2000) on page 81.
%
% example:
% f = plateRectModes(3, 3, 1, 1);

% total number of modes in each direction (number of half waves)
% this is just what makes sense for now. i might change this depending on
% resulting sounds
Nm = floor(dx/2) + 1;
Nn = floor(dy/2) + 1;

% store the modal frequencies in f
f = zeros(Nm, Nn);

for m=1:Nm
    for n=1:Nn
        f(m, n) = 0.453 * cL * h * (((m+2)/dx)^2 + ((n+2)/dy)^2);
    end
end

% normalize frequencies
f = f/f(1);

end


