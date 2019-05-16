% pitchGlideExamples.m
% This script shows us how to:
%   1. set an exponential decaying pitch glide
%   2. control how quickly the pitch glide decays
%   3. set up a linearly increasing pitch glide? or squareroot?

addpath(genpath('../proofOfConcept'));

savePlots = 1;
plotOutDir = 'figures/';

% input parameters
fs = 44100;
dur = 0.9;
decayT60 = 0.6;     % time i want the pitch glide envelope to be at -60dB
B = 0.5;           % timbre control

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

%w0 = 2*pi*f0_0;


%% Linear increase

f0_0 = 140;         % starting frequency, Hz
f0_1 = 840;          % ending frequency, Hz

m = (f0_1 - f0_0)/dur;
b = f0_0;
f0Lin = m*nT + b;
w0Tilde = 2*pi*f0Lin;

%% LOOPBACK FM LINEAR INCREASE

wc = w0Tilde(end)/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0Tilde/wc).^2);

yLin = zeros(1, N);
yLin(1) = 1;
for i=2:N
    yLin(i) = exp(j*wc*T*(1 + BTilde(i) * real(yLin(i-1)))) * yLin(i-1);
end

%% CLOSED FORM LOOPBACK FM LINEAR INCREASE

% Theta0 = integral[2*pi * (m*nT + b) dnT]
Theta0 = 2*pi * (m/2 * nT.^2 + b*nT);
b0 = (sqrt(1 - B^2) - 1)/B;
YLin = (b0 + exp(1j.*Theta0)) ./ (1 + b0.*exp(1j.*Theta0));


%% linearly increasing pitch glide plots

figure

subplot(211)
spectrogram(real(yLin), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(n*T*1000, f0Lin/1000, 'r', 'linewidth', 2)
ylim([0 1])
title('Pitch glide increasing linearly with z_c(n)');
colorbar('off')

subplot(212)
spectrogram(real(YLin), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(n*T*1000, f0Lin/1000, 'r', 'linewidth', 2)
ylim([0 1])
title('Pitch glide increasing linearly with z_0(n)');
colorbar('off')

if savePlots
    H = figure
    spectrogram(real(YLin), hann(1024), 512, 2048, fs, 'yaxis')
    hold on
    plot(n*T*1000, f0Lin/1000, 'r', 'linewidth', 2)
    ylim([0 1])
    title('Pitch glide increasing linearly with z_0(n)');
    set(gca, 'fontsize', 15)
    fig = H
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.6];
    print([plotOutDir 'pitchGlideLinIncrease'], '-depsc', '-r0')
    %close
end



%% EXPONENTIAL DECAY


f0_0 = 840;         % starting frequency, Hz
f0_1 = 140;          % ending frequency, Hz

% exponentially decaying envelope according to the T60
A = 1;
tau = -(decayT60)/log(0.001);
e = A * exp(-nT./tau);

% scaled f0 function
f0Exp = (f0_0 - f0_1) * e + f0_1;

figure
plot(nT, f0Exp);
hold on
plot([decayT60 decayT60], [f0_0 f0_1], 'k--')
plot([0 decayT60], [f0_1 f0_1], 'k--')


%% LOOPBACK FM EXPONENTIAL DECAY

% loopback FM variables
% we need to make sure the w0 <= wc. w0Tilde(1) will always be the largest
% value since this is an exponentially decreasing function. We can use
% w0Tilde(1) and a user-given B to solve for wc.
w0Tilde = 2*pi*f0Exp;
wc = w0Tilde(1)/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0Tilde/wc).^2);

yExp = zeros(1, N);
yExp(1) = 1;
for i=2:N
    yExp(i) = exp(j*wc*T*(1 + BTilde(i) * real(yExp(i-1)))) * yExp(i-1);
end

%soundsc(real(yDec) .* env, fs)

%% CLOSED FORM LOOPBACK FM EXPONENTIAL DECAY

% Theta0 = integral [(f0_0 - f0_1) * A * exp(-nT/tau) + f0_1 dnT] <-- that's the equation
Theta0 = (2 * pi * (f0_0 - f0_1)) * (-tau * A) * exp(-nT/tau) + (2 * pi * f0_1 * nT); 

b0 = (sqrt(1 - B^2) - 1)/B;
YExp = (b0 + exp(1j.*Theta0)) ./ (1 + b0.*exp(1j.*Theta0));

%% decaying pitch glide plots

figure

subplot(211)
spectrogram(real(yExp), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(n*T*1000, f0Exp/1000, 'r', 'linewidth', 2)
ylim([0 1])
title('Pitch glide decreasing exponentially with z_c(n)');

subplot(212)
spectrogram(real(YExp), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(n*T*1000, f0Exp/1000, 'r', 'linewidth', 2)
ylim([0 1])
title('Pitch glide decreasing exponentially with z_0(n)');

if savePlots
    H = figure
    spectrogram(real(YExp), hann(1024), 512, 2048, fs, 'yaxis')
    hold on
    plot(n*T*1000, f0Exp/1000, 'r', 'linewidth', 2)
    ylim([0 1])
    colorbar('off')
    title('Pitch glide decreasing exponentially with z_0(n)');
    set(gca, 'fontsize', 15)
    fig = H
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.6];
    print([plotOutDir 'pitchGlideExpDecay'], '-depsc', '-r0')
    %close
end

%% SQUARE ROOT INCREASE

f0_0 = 140;
f0_1 = 840;

% square-root increase type envelope
c = f0_0;
%a = (f0_1 - f0_0)/sqrt(T60Samp);
%f0Tilde = a*sqrt(n) + c;
a = (f0_1 - f0_0)/sqrt(decayT60);
f0Sqrt = a*sqrt(nT) + c;

figure
plot(nT, f0Sqrt);
hold on
plot([decayT60 decayT60], [f0_0 f0_1], 'k--')
plot([0 decayT60], [f0_1 f0_1], 'k--')


%% LOOPBACK FM SQUAREROOT INCREASE

% loopback FM variables
% we need to make sure the w0 <= wc. w0Tilde(N-1) will always be the largest
% value since this is an increasing square-root function. We can use
% w0Tilde(N-1) and a user-given B to solve for wc.
w0Tilde = 2*pi*f0Sqrt;
wc = w0Tilde(N-1)/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0Tilde/wc).^2);

ySqrt = zeros(1, N);
ySqrt(1) = 1;
for i=2:N
    ySqrt(i) = exp(j*wc*T*(1 + BTilde(i) * real(ySqrt(i-1)))) * ySqrt(i-1);
end

%% STRETCHED APF SQUAREROOT INCREASE

% Theta0 = integral[2*pi * (a*sqrt(nT) + c) dnT]
Theta0 = 2*pi * ((2/3) * a * (nT.^(3/2)) + c*nT);
b0 = (sqrt(1 - B^2) - 1)/B;
YSqrt = (b0 + exp(1j.*Theta0)) ./ (1 + b0.*exp(1j.*Theta0));


%% SQUAREROOT INCREASING PITCH GLIDE

figure

subplot(211)
spectrogram(real(ySqrt), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(n*T*1000, f0Sqrt/1000, 'r', 'linewidth', 2)
ylim([0 1])
title('Pitch glide increasing by a square root function with z_c(n)');
colorbar('off')

subplot(212)
spectrogram(real(YSqrt), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(n*T*1000, f0Sqrt/1000, 'r', 'linewidth', 2)
ylim([0 1])
title('Pitch glide increasing by a square root function with z_0(n)');
colorbar('off')

if savePlots
    H = figure
    spectrogram(real(YSqrt), hann(1024), 512, 2048, fs, 'yaxis')
    hold on
    plot(n*T*1000, f0Sqrt/1000, 'r', 'linewidth', 2)
    ylim([0 1])
    colorbar('off')
    title('Pitch glide increasing by a square root function with z_0(n)');
    set(gca, 'fontsize', 15)
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.6];
    print([plotOutDir 'pitchGlideSqrtIncrease'], '-depsc', '-r0')
    %close
end

%% with the setup we already have, how does w0 change with different B
% functions?

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



