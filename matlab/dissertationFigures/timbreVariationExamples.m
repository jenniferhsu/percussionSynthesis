addpath(genpath('../proofOfConcept'));

savePlots = 1;
plotOutDir = 'figures/';

% input parameters
fs = 44100;
dur = 0.9;
decayT60 = 0.6;     % time i want the pitch glide envelope to be at -60dB
B = 0.5;            % timbre control
f0 = 1000;

% derived parameters
N = dur*fs;
T = 1/fs;
n = 0:1:N-1;
nT = 0:T:(dur-T);
w0 = 2*pi*f0;

% decay envelope
A_e = 1;
tau_e = -(N-1)/log(0.001);
env = A_e * exp(-n/tau_e);


%% CLOSED FORM LOOPBACK FM LINEARLY INCREASE b0

b00 = 0;
b01 = 0.8;
k = (b01 - b00)/N;
l = 1 - k*N;
b0 = k * n + l;
y0Lin = (b0 + exp(1j*w0*n*T)) ./ (1 + b0.*exp(1j*w0*n*T));


%% LOOPBACK FM LINEAR INCREASE

BTilde = - 2*b0 ./(b0.^2 + 1);
ycLin = zeros(1, N);
ycLin(1) = 1;
for i=2:N
    ycLin(i) = exp(1j*wc*T*(1 + BTilde(i) * real(ycLin(i-1)))) * ycLin(i-1);
end



%% linearly increasing pitch glide plots

figure

subplot(211)
spectrogram(real(y0Lin), hann(1024), 512, 2048, fs, 'yaxis')
%ylim([0 1])
xlim([1 dur*100])
title('Timbre increasing linearly with z_0(n)');
colorbar('off')

subplot(212)
spectrogram(real(ycLin), hann(1024), 512, 2048, fs, 'yaxis')
%ylim([0 1])
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
