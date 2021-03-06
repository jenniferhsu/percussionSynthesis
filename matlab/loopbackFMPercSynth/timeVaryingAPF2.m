function yApf = timeVaryingAPF2(x, f_pi, f_b, M, f_m, fs)
%%timeVaryingAPF2(x, f_pi, f_b, M, f_m, fs)
% Implements the time-varying 2nd order all-pass filter from Greg Surges's
% paper "Generative Feedback Networks using Time-Varying Allpass Filters"
%
% inputs:
%   x: the signal to be filtered
%   f_pi: the frequency at which the phase is -pi
%       this is also the center frequency
%   f_b: the bandwidth of the generated sidebands
%   M: modulation depth
%   f_m: frequency of modulation
%   fs: sampling rate
%
% output:
%   yApf: the output signal of the same size as input x
%
% author: Jennifer Hsu
% date: Spring 2019

N = length(x);
nT = (0:(N-1))/fs;
yApf = zeros(size(x));

% time-varying intermediate parameters
f_pi_tilde = f_pi + M*cos(2*pi*f_m*nT);
d_tilde = -cos(2*pi*f_pi_tilde/fs);

% filter memory
z = zeros(1,2); 

% filter the signal sample-by-sample
for n=1:N
    
    d = d_tilde(n);
    c = (tan((pi*f_b)/fs) - 1)/(tan((pi*f_b)/fs) + 1);
    
    r1 = acos(-c);
    r2 = acos(-d);
    
    % matrix method
    zOld = z;
    yApf(n) = cos(r1)*x(n) - sin(r1)*cos(r2)*zOld(1) + sin(r1)*sin(r2)*zOld(2);
    z(1) = sin(r1)*x(n) + cos(r1)*cos(r2)*zOld(1) - sin(r2)*cos(r1)*zOld(2);
    z(2) = sin(r2)*zOld(1) + cos(r2)*zOld(2);

end