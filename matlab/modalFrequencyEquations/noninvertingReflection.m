% this script calculates the modal frequencies of a rectangular
% plate with "simply-supported" edges -- I think this means that we
% have non-inverting reflections on all four sides.
%
% the equation for the modal frequencies is from The Science of
% Percussion Instruments (2000) by Tom Rossing on page 81

% load basic time-domain functions to compare with
addpath(genpath('../basicMeshFunctions/'));

% mesh size
Nj = 3;
dx = Nj;
dy = Nj;

h = 1;	% plate thickness
cL = 1;	% speed of longitudinal waves

% mode numbers - number of half waves 
Nm = floor(dx/2)+1;
Nn = floor(dy/2)+1;

% modal frequencies
f = zeros(Nm, Nn);

for m=1:Nm
    for n=1:Nn
        % it looks like i need to add an extra 1 to the nodal numbers?
    	f(m,n) = 0.453 * cL * h * (((m+2)/dx)^2 + ((n+2)/dy)^2);
    end
end

% normalize the frequencies
f = f/f(1);

% compare with basic mesh time-domain implementation
fs = 44100;
refl = 1;
N = fs;
y = mesh2DReflCoeff(1, Nj, N, fs, Nj, Nj, Nj, Nj, refl, 0);

% get the FFT of time-domain implementation
Nfft = 1024;
fAxis = (fs/2)*linspace(0, 1, (Nfft/2)+1);

Y = fft(y, Nfft);
YAbs = abs(Y);
YAbsPos = YAbs(1:(Nfft/2)+1);

% get the fundamental frequency of the mesh IR
[pks, locs] = findpeaks(YAbsPos);
f0 = fAxis(locs(1));

% see if our calculated modal frequencies match up
freqPeaks = unique(f(:)*f0)

figure
plot(fAxis, YAbsPos)
hold on
for i=1:length(freqPeaks)
    plot([freqPeaks(i) freqPeaks(i)], [0 max(YAbsPos)], 'r--');
 end
ylim([0 max(YAbsPos)])

