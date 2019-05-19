% Chapter5Plots.m
% This script plots all the figures from Chapter 5.
%
% author: Jennifer Hsu
% date: Spring 2019

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
T60 = 0.8;

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
kick0c = loopbackFMzc(fc, B, BEnd, g, BGlideType, dur, fs);
kick0c = kick0c .* w;

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
sgtitle('Kick drum example: static pitch and timbre')
if savePlots
    saveas(gcf, [figDir filePrefix 'kickStaticPitchAndTimbre'], 'epsc')
end

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_kick_staticPitchAndTimbre' '.wav'], ...
                scaleForSavingAudio(real(kick0)), fs);
end

%% Figure 5.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Exponentially decreasing pitch glide: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 1.0;
T = 1/fs;
N = dur*fs;
n = 0:N-1;

% exponentially decreasing pitch glide parameters
fx = 100;
fy = 40;
r = 0.001;
A = 1;
t_d = 0.6;

n_d = t_d * fs;

% calculate the exponentially decreasing pitch glide from 0 to 1
tau = -(n_d*T)/log(r/A);
d = exp(-n*T/tau);

% scale using fx and fy
f0n = (fx - fy) * d + fy;

figure
plot(n*T, f0n, 'linewidth', 2)
xlabel('Time (seconds)');
ylabel('Amplitude (linear)');
title('Kick drum example: exponentially decreasing pitch glide f_{0,0}(n)');
grid on
set(gca, 'FontSize', 15);
if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, [figDir filePrefix 'KickDrumf0Tilde'], 'epsc')
end


%% Figure 5.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Exponentially decreasing pitch glide B(n) for z_{c,i}(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0n = 2*pi*f0n;
wc = w0n(1);
Bn = sqrt(1 - (w0n./wc).^2);

figure
plot(n*T, Bn, 'linewidth', 2)
xlabel('Time (seconds)');
ylabel('B(n)');
title('Kick drum example: exponentially decreasing pitch glide B(n) function');
grid on
set(gca, 'FontSize', 15);
if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, [figDir filePrefix 'KickDrumBn'], 'epsc')
end


%% Figure 5.4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Exponentially decreasing pitch glide kick drum synthesis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loopback FM signal (z0)
fs = 44100;
dur = 1;
f0 = 100;
f0End = 40;
b0 = 0.2;
pitchGlideType = 'exp';
aStruct = struct();
aStruct.T60 = 0.6;
kick1 = loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, aStruct);

% envelope
A0 = 1;
T60 = 0.8;

T = 1/fs;
N = dur*fs;
n = 0:N-1;

n60 = T60*fs;
A60 = A0/10^(60/20);
tau = -n60*T/log(A60/A0);

w = A0*exp(-n*T/tau);
kick1 = kick1 .* w;

% loopback FM (zc)
w0n = 2*pi*f0n;
wc = w0n(1);
Bn = sqrt(1 - (w0n./wc).^2);

kick1c = loopbackFMzc(fc, Bn, 0, 0, 'useB', dur, fs);
kick1c = kick1c .* w;

% plot
figure
subplot(211)
spectrogram(real(kick1), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(n*T*1000, f0n/1000, 'r', 'linewidth', 2)
colorbar('off');
title('z_{0,0}(n) oscillator');
xlim([0 1000])
ylim([0 1])
set(gca, 'FontSize', 15);
set(gca, 'XTick', [200 400 600 800 1000], 'XTickLabel', [200 400 600 800 1000])

subplot(212)
spectrogram(real(kick1c), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(n*T*1000, f0n/1000, 'r', 'linewidth', 2)
colorbar('off');
title('z_{c,0}(n) oscillator');
xlim([0 1000])
ylim([0 1])
set(gca, 'FontSize', 15);
set(gca, 'XTick', [200 400 600 800 1000], 'XTickLabel', [200 400 600 800 1000])
sgtitle('Kick drum example: exponentially decreasing pitch glide')
if savePlots
    saveas(gcf, [figDir filePrefix 'kickPitchGlide'], 'epsc')
end

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_kick_pitchGlidez0' '.wav'], ...
                scaleForSavingAudio(real(kick1)), fs);
    audiowrite([audioDir 'hsu_kick_pitchGlidezc' '.wav'], ...
                scaleForSavingAudio(real(kick1c)), fs);
end

%% Figure 5.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Exponentially decreasing pitch glide kick drum synthesis
% and linearly decreasing timbre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loopback FM signal (z0)
fs = 44100;
dur = 1;
f0 = 100;
f0End = 40;
N = fs*dur;
b0 = linspace(0.5, 0, N);
pitchGlideType = 'exp';
aStruct = struct();
aStruct.T60 = 0.6;
kick2 = loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, aStruct);

% envelope
A0 = 1;
T60 = 0.8;

T = 1/fs;
N = dur*fs;
n = 0:N-1;

n60 = T60*fs;
A60 = A0/10^(60/20);
tau = -n60*T/log(A60/A0);

w = A0*exp(-n*T/tau);
kick2 = kick2 .* w;

% plot
figure
spectrogram(real(kick2), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(n*T*1000, f0n/1000, 'r', 'linewidth', 2)
colorbar('off');
title('z_{0,0}(n) oscillator with pitch glide and timbre variation');
xlim([0 1000])
ylim([0 1])
set(gca, 'FontSize', 15);
set(gca, 'XTick', [200 400 600 800 1000], 'XTickLabel', [200 400 600 800 1000])

if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, [figDir filePrefix 'kickPitchGlideTimbreVariation'], 'epsc')
end

if saveAudio
    audiowrite([audioDir 'hsu_kick_pitchAndTimbreVariation' '.wav'], ...
                scaleForSavingAudio(real(kick2)), fs);
end

%% Figure 5.6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Applying time-varying allpass filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loopback FM signal (z0)
fs = 44100;
dur = 1;
f0 = 100;
f0End = 40;
N = fs*dur;
b0 = linspace(0.5, 0, N);
pitchGlideType = 'exp';
aStruct = struct();
aStruct.T60 = 0.6;
kick3 = loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, aStruct);

% envelope
A0 = 1;
T60 = 0.8;

T = 1/fs;
N = dur*fs;
n = 0:N-1;

n60 = T60*fs;
A60 = A0/10^(60/20);
tau = -n60*T/log(A60/A0);

w = A0*exp(-n*T/tau);

% time-varying allpass filter parameters
TVAPFParams.fmVec = [100];
TVAPFParams.MVec = [120];
TVAPFParams.fbVec = [2500];
TVAPFParams.fpiVec = [f0];

[kick3, yMat] = applyTimeVaryingAPF2(kick3, w, fs, TVAPFParams);

% plot
figure
spectrogram(real(kick3), hann(256), 128, 1024, fs, 'yaxis');
%hold on
%plot(n*T*1000, f0n/1000, 'r', 'linewidth', 2)
colorbar('off');
title('z_{0,0}(n) kick drum with time-varying allpass filter');
xlim([0 1000])
ylim([0 1])
set(gca, 'FontSize', 15);
set(gca, 'XTick', [200 400 600 800 1000], 'XTickLabel', [200 400 600 800 1000])

if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, [figDir filePrefix 'kickTimeVaryingAPF'], 'epsc')
end

if saveAudio
    audiowrite([audioDir 'hsu_kick_timeVaryingAPF' '.wav'], ...
                scaleForSavingAudio(real(kick3)), fs);
end

%% Figure 5.7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Applying Commuted Synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resIRWav = '../loopbackFMPercSynth/resonatorIRs/taiko/taiko2_fixedOffset.wav';
excitationType = 'rc';
excitationParams = struct();
excitationParams.winLength = 4;
kick4 = applyCommutedSynthesis(kick3, resIRWav, excitationType, excitationParams, fs);

% plot
figure
spectrogram(real(kick4), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('z_{0,0}(n) kick drum with commuted synthesis');
ylim([0 1])
set(gca, 'FontSize', 15);
if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, [figDir filePrefix 'kickCommutedSynthesis'], 'epsc')
end

if saveAudio
    audiowrite([audioDir 'hsu_kick_commutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(kick3)), fs);
end

%% Figure 5.8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snare Drum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 1.0;

% === snare - set up argument struct ===
argStruct = struct();
argStruct.sinusoidArgs = struct();
argStruct.zcArgs = struct();
argStruct.z0Args = struct();

% === snare - set up envelope ===
e0Vec = [1, 0.2]';
T60Vec = [0.6, 0.5]';
env = envMat(e0Vec, T60Vec, dur, fs);
wg = env(1,:)';  % global amplitude envelope

Nf = length(e0Vec);
N = dur*fs;

% === snare synthesis - loopback FM parameters ===

% traditional MS parameters
argStruct.sinusoidArgs.f0Vec = [185, 330];
argStruct.sinusoidArgs.f0EndVec = [185, 330];
argStruct.sinusoidArgs.pitchGlideTypeVec = {'none', 'none'};
argStruct.sinusoidArgs.zcArgsVec = [0, 0];

% loopback FM z0 parameters
argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = [185, 330];
argStruct.z0Args.f0EndVec = [140, 280];
argStruct.z0Args.pitchGlideTypeVec = {'exp', 'exp'};
for f=1:Nf
    argStruct.z0Args.b0Mat(f,:) = linspace(0.999*((Nf-(f-1))/Nf), 0.001, N);
    argStruct.z0Args.b0Mat(f,:) = ones(1,N) * 0.5;
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f) = struct();
end
argStruct.z0Args.zcArgsVec(1).T60 = 0.6;
argStruct.z0Args.zcArgsVec(2).T60 = 0.6;

% === snare synthesis ===

% snare - traditional MS
[snare_s, snare_s_Mat] = loopbackFMMS('s', env, argStruct, fs);

% snare - loopback FM z0
[snare_z0, snare_z0_Mat] = loopbackFMMS('z0', env, argStruct, fs);

% === time-varying allpass filter ===

TVAPFParams.fpiVec = argStruct.sinusoidArgs.f0Vec;
TVAPFParams.fbVec = [100 200];
TVAPFParams.MVec = [1000 250];
TVAPFParams.fmVec = [300 2000];

% time-varying APF - traditional MS
[ySnare_s, ySnare_s_Mat] = applyTimeVaryingAPF2(snare_s_Mat, env, fs, TVAPFParams);

% time-varying APF - loopback FM z0
[ySnare_z0, ySnare_z0_Mat] = applyTimeVaryingAPF2(snare_z0_Mat, env, fs, TVAPFParams);

% === commuted synthesis ===

% commuted synthesis parameters
resIRWav = '../loopbackFMPercSynth/resonatorIRs/impulse.wav';
excitationType = 'nb';
excitationParams = struct();
excitationParams.durNB = 0.6;
excitationParams.lowFreq = 120;
excitationParams.highFreq = 8000;

% snare - traditional MS commuted synthesis
snare_s_CS = applyCommutedSynthesis(real(snare_s), resIRWav, excitationType, excitationParams, fs);
snare_s_CS = snare_s_CS .* wg;

% snare - loopback FM z0 commuted synthesis
snare_z0_CS = applyCommutedSynthesis(real(snare_z0), resIRWav, excitationType, excitationParams, fs);
snare_z0_CS = snare_z0_CS .* wg;

% snare - traditional MS with time-varying APF - commuted synthesis
ySnare_s_CS = applyCommutedSynthesis(real(ySnare_s), resIRWav, excitationType, excitationParams, fs);
ySnare_s_CS = ySnare_s_CS .* wg;

% snare - loopback FM z0 with time-varying APF - commuted synthesis
ySnare_z0_CS = applyCommutedSynthesis(real(ySnare_z0), resIRWav, excitationType, excitationParams, fs);
ySnare_z0_CS = ySnare_z0_CS .* wg;

% === plot ===
figure
subplot(411)
spectrogram(real(snare_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS');
xlim([0 600])
ylim([0 10])
set(gca, 'FontSize', 15);

subplot(412)
spectrogram(real(snare_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS');
xlim([0 600])
ylim([0 10])
set(gca, 'FontSize', 15);

subplot(413)
spectrogram(real(ySnare_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS with time-varying allpass filters');
xlim([0 600])
ylim([0 10])
set(gca, 'FontSize', 15);

subplot(414)
spectrogram(real(ySnare_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with time-varying allpass filters');
xlim([0 600])
ylim([0 10])
set(gca, 'FontSize', 15);

sgtitle('Snare drum synthesis example', 'FontSize', 15)

if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    saveas(gcf, [figDir filePrefix 'snare'], 'epsc')
end

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_snare_traditionalMS' '.wav'], ...
                scaleForSavingAudio(real(snare_s_CS)), fs);
    audiowrite([audioDir 'hsu_snare_loopbackFM' '.wav'], ...
                scaleForSavingAudio(real(snare_z0_CS)), fs);
    audiowrite([audioDir 'hsu_snare_traditionalMSTimeVaryingAPF' '.wav'], ...
                scaleForSavingAudio(real(ySnare_s_CS)), fs);
    audiowrite([audioDir 'hsu_snare_loopbackFMTimeVaryingAPF' '.wav'], ...
                scaleForSavingAudio(real(ySnare_z0_CS)), fs);
end