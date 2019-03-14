% pitchGlideExamples.m
% This script shows us how to:
%   1. set an exponential decaying pitch glide
%   2. control how quickly the pitch glide decays
%   3. set up a linearly increasing pitch glide? or squareroot?

addpath(genpath('../proofOfConcept'));

% input parameters
fs = 44100;
dur = 0.9;
decayT60 = 0.5;     % time i want the pitch glide envelope to be at -60dB
f0_0 = 100;         % starting frequency, Hz
f0_1 = 40;          % ending frequency, Hz
B = 0.01;           % timbre control

% derived parameters
N = dur*fs;
T = 1/fs;
n = 0:1:N-1;
nT = 0:T:(dur-T);
T60Samp = decayT60*fs;

% decay envelope
A_e = 1;
tau_e = -(N-1)/log(0.001);
env = A_e * exp(-n/tau_e);

w0 = 2*pi*f0_0;

%% loopback FM - exponential decay 1

% exponentially decaying envelope according to the T60
A = 1;
tau = -(T60Samp)/log(0.001);
e = A * exp(-n./tau);

% scaled f0 function
f0Tilde = (f0_0 - f0_1) * e + f0_1;

% loopback FM variables
% we need to make sure the w0 <= wc. w0Tilde(1) will always be the largest
% value since this is an exponentially decreasing function. We can use
% w0Tilde(1) and a user-given B to solve for wc.
w0Tilde = 2*pi*f0Tilde;
wc = w0Tilde(1)/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0Tilde/wc).^2);

yExp = zeros(1, N);
yExp(1) = 1;
for i=2:N
    yExp(i) = exp(j*wc*T*(1 + BTilde(i) * real(yExp(i-1)))) * yExp(i-1);
end

%spectrogram(real(yExp), hann(256), 128, 1024, fs, 'yaxis')
spectrogram(real(yExp), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(n*T*1000, f0Tilde/1000, 'r')
ylim([0 1])
title('loopback FM - exponential w0 pitch glide');


%soundsc(real(yExp) .* env, fs)
%audiowrite('audioExamples/kickExp.wav', scaleForSavingAudio(real(yExp) .* env), fs);

%% square root-type increase

wc = 2*pi*220;

% linear increase of B is a squareroot decrease of w0
BB1 = 0:0.1:1;
x1 = wc * sqrt(1 - BB1.^2);

% linear decrease of B is a squareroot increase of w0
BB2 = 1:-0.1:0;
x2 = wc * sqrt(1 - BB2.^2);

% exponential decrease of B 
A = 1;
tau = -(0.5*fs)/log(0.001);
e = A * exp(-n./tau);
BB3 = e;
x3 = wc * sqrt(1 - BB3.^2);

A = 1;
tau = -(0.25*fs)/log(0.001);
e = A * exp(-n./tau);
BB4 = e;
x4 = wc * sqrt(1 - BB4.^2);

plot(x3)
hold on
plot(x4)
% cool, so the steeper the decrease of B, the steeper the increase of w0

% for increasing pitch glides, make sure that wc is equal to or larger than
% the final/last value of w0 for the pitch glide.


