% Chapter5Plots.m
% This script plots all the figures from Chapter 5.

addpath(genpath('../loopbackFMPercSynth/'));
addpath(genpath('../helperFunctions/'));
savePlots = 1;
saveAudio = 1;

figDir = 'figures/';
audioDir = 'soundExamples/';
filePrefix = 'Chapter5_';

%% Figure 5.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Static Pitch and Timbre: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loopback FM signal (z0)
fs = 44100;
dur = 1.0;
f0 = 100;
f0End = 100;
b0 = 0.2;
pitchGlideType = 'none';
kick0 = loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs);

% envelope
A0 = 1;
T60 = 0.6;

T = 1/fs;
N = dur*fs;
n = 0:N-1;

n60 = T60*fs;
A60 = A0/10^(60/20);
tau = -n60*T/log(A60/A0);

w = A0*exp(-n*T/tau);
kick0 = kick0 .* w;

% loopback FM (zc)
B = -2*b0/(b0^2 + 1);
BEnd = B;
BGlideType = 'none';
g = 0;
fc = f0/sqrt(1 - B^2);
kick01 = loopbackFMzc(fc, B, BEnd, g, BGlideType, dur, fs);
kick01 = kick01 .* w;

% plot
figure
subplot(211)
plot(n*T*1000, real(kick0));
xlabel('Time (ms)');
ylabel('Amplitude (linear)');
title('Time-domain signal');
set(gca, 'FontSize', 15);
grid on
subplot(212)
spectrogram(real(kick0), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Spectrogram');
xlim([0 1000])
ylim([0 2])
set(gca, 'FontSize', 15);
set(gca, 'XTick', [200 400 600 800 1000], 'XTickLabel', [200 400 600 800 1000])
sgtitle('Kick drum synthesis with static pitch and timbre')
if savePlots
    saveas(gcf, [figDir filePrefix 'kickStaticPitchAndTimbre'], 'epsc')
end

% save audio
if saveAudio
    audiowrite([audioDir filePrefix 'kickStaticPitchAndTimbre' '.wav'], ...
                scaleForSavingAudio(real(kick0)), fs);
end