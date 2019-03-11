%% This script shows how we can do an exponentially and linearly 
% decaying pitch glide by setting w0 (and not B!). The sound here is really
% similar to a kick drum.

% input parameters
fs = 44100;
dur = 0.6;
alpha = 0.00025;
f0_1 = 100;         % starting frequency, Hz
f0_2 = 40;          % ending frequency, Hz
B = 0.01;           % timbre control

% derived parameters
N = dur*fs;
T = 1/fs;
n = 0:1:N-1;
nT = 0:T:(dur-T);

% decay envelope
A_e = 1;
tau_e = -(N-1)/log(0.001);
env = A_e * exp(-n/tau_e);

w0 = 2*pi*f0_1;


%% loopback FM - exponential

tau = -(N-1)/log(f0_2/f0_1);
w0Tilde = w0 * exp(-n./tau);
wc = w0/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0Tilde/wc).^2);

yExp = zeros(1, N);
yExp(1) = 1;
for i=2:N
    yExp(i) = exp(j*wc*T*(1 + BTilde(i) * real(yExp(i-1)))) * yExp(i-1);
end

spectrogram(real(yExp), hann(256), 128, 1024, fs, 'yaxis')
ylim([0 2])
title('loopback FM - exponential w0 pitch glide');

%soundsc(real(yExp) .* env, fs)

%% loopback FM - linear

m = 2 * pi * (-60/dur);
b = 2 * pi * 100;
w0Tilde = m*nT + b;
wc = w0/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0Tilde/wc).^2);

yLin = zeros(1, N);
yLin(1) = 1;
for i=2:N
    yLin(i) = exp(j*wc*T*(1 + BTilde(i) * real(yLin(i-1)))) * yLin(i-1);
end

spectrogram(real(yLin), hann(256), 128, 1024, fs, 'yaxis')
ylim([0 2])
title('loopback FM - linear w0 pitch glide');

%soundsc(real(yLin) .* env, fs)


%% stretched APF - exponential

% exponential decrease for w0
A = w0;
tau = -dur/log(2*pi*40/w0);
b0 = (sqrt(1 - B^2) - 1)/B;
%ThetaH = integral [exp(-nT/tau dn] <-- that's the equation
ThetaH = (-tau * A) * exp(-nT/tau); 
YExp = (b0 + exp(1j.*ThetaH)) ./ (1 + b0.*exp(1j.*ThetaH));

spectrogram(real(YExp), hann(256), 128, 1024, fs, 'yaxis')
ylim([0 2])
title('stretched APF - exponential w0 pitch glide');


%% stretched APF - linear 
m = 2 * pi * (-60/dur);
b = 2 * pi * 100;
b0 = (sqrt(1 - B^2) - 1)/B;
%ThetaH = integral [m * nT + b d(nT)] <-- that's the equation (why is T
%included here, but not in the version above?
ThetaH = (m*nT.^2)/2 + b*nT;
YLin = (b0 + exp(1j.*ThetaH)) ./ (1 + b0.*exp(1j.*ThetaH));

spectrogram(real(YLin), hann(256), 128, 1024, fs, 'yaxis')
ylim([0 2])
title('stretched APF - linear w0 pitch glide');