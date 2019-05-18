function env = envMat(e0Vec, T60Vec, dur, fs)
%envMat(e0Vec, T60Vec) creates an amplitude envelope matrix to be used with 
%   loopback FM MS.  env has size [Nf, N] where      
%       Nf = number of modal frequencies
%       N = desired length of signal in samples
%
% inputs:
%   e0Vec: vector of initial amplitude values of size [Nf, 1]
%   T60Vec: vector of T60 times in seconds of size [Nf, 1]
%   dur: desired duration of signal in seconds
%   fs: sample rate (samples/sec)
%   
% output:
%   env: envelope matrix of size [Nf,N] where
%
% example:
% e0Vec = [1, 0.9, 0.89, 0.8, 0.75]';
% T60Vec = [0.9, 0.88, 0.87, 0.85, 0.77]';
% env = envMat(e0Vec, T60Vec, dur, fs);

N = dur*fs;
Nf = length(e0Vec);
T = 1/fs;

n60Vec = T60Vec * fs;           % T60 in samples
e60Vec = e0Vec ./ 10^(60/20);   % amplitude envelope end value

A = e0Vec;
tau = (-n60Vec * T) ./ log(e60Vec ./ A);

nT = linspace(0, dur, N);
nTMat = repmat(nT, Nf, 1);

env = A .* exp(-nTMat ./ tau);


end