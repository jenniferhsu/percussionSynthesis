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
    %argStruct.z0Args.b0Mat(f,:) = linspace(0.999*((Nf-(f-1))/Nf), 0.001, N);
    argStruct.z0Args.b0Mat(f,:) = ones(1,N) * 0.5;
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f) = struct();
end
argStruct.z0Args.zcArgsVec(1).T60 = 0.6;
argStruct.z0Args.zcArgsVec(2).T60 = 0.6;

% === snare synthesis - traditional MS and loopback FM MS ===

% snare - traditional MS
[snare_s, snare_s_Mat] = loopbackFMMS('s', env, argStruct, fs);

% snare - loopback FM z0
[snare_z0, snare_z0_Mat] = loopbackFMMS('z0', env, argStruct, fs);

% === snare synthesis - time-varying allpass filter ===

TVAPFParams.fpiVec = argStruct.sinusoidArgs.f0Vec;
TVAPFParams.fbVec = [100 200];
TVAPFParams.MVec = [1000 250];
TVAPFParams.fmVec = [300 2000];

% time-varying APF - traditional MS
[ySnare_s, ySnare_s_Mat] = applyTimeVaryingAPF2(snare_s_Mat, env, fs, TVAPFParams);

% time-varying APF - loopback FM z0
[ySnare_z0, ySnare_z0_Mat] = applyTimeVaryingAPF2(snare_z0_Mat, env, fs, TVAPFParams);

% === snare synthesis - commuted synthesis ===

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

% === snare synthesis - plot ===
figure
subplot(411)
spectrogram(real(snare_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS with CS');
xlim([0 600])
ylim([0 5])
set(gca, 'FontSize', 15);

subplot(412)
spectrogram(real(snare_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with CS');
xlim([0 600])
ylim([0 5])
set(gca, 'FontSize', 15);

subplot(413)
spectrogram(real(ySnare_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS with CS and time-varying APFs');
xlim([0 600])
ylim([0 5])
set(gca, 'FontSize', 15);

subplot(414)
spectrogram(real(ySnare_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS w/ CS and time-varying APFs');
xlim([0 600])
ylim([0 5])
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



%% Figure 5.9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Marimba
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 2;

% === derived parameters === 
% feedback FM
B = 0.3;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1);             % feedback FM pitch glide coefficient

% === ideal bar modal frequencies === 
% from Science of Percussion Instruments page 47

f1 = 440;
Nf = 7;
fVecBar = zeros(1, Nf);
for i=1:Nf
    if i==1
        fVecBar(i) = 3.011^2;
    else
        fVecBar(i) = (2*i+1)^2;
    end
end
fVecBar = fVecBar/fVecBar(1);
fVecBar = f1 * fVecBar;

% center frequencies for feedback FM if sounding frequencies = membrane
% modal freqeuncies
fcVecBar = fVecBar./(sqrt(1-B^2));

% === set up decay envelopes ===
% want the lowest modal frequencies to have a larger starting amplitude 
% than the highest modal frequencies
aStart = zeros(1, Nf);
aStart(1) = 1;              % highest amplitude for lowest modal freq
aStart(Nf) = 0.01;          % lowest amplitude for highest modal freq
tau = -(Nf-1) / log(aStart(Nf));
a0 = 1/exp(-1/tau);
aStart(2:Nf-1) = a0*exp(-(2:Nf-1)/tau);     % exponentially decreasing 
                                            % starting amplitudes

e = g.^(linspace(0, N, N));   % exponential decay

env = zeros(Nf, N);
for i=1:Nf
    env(i,:) = aStart(i) * e;
end
wg = env(1,:)';

% === marimba synthesis - loopback FM parameters ===

% traditional MS parameters
argStruct.sinusoidArgs.f0Vec = fVecBar;
argStruct.sinusoidArgs.f0EndVec = fVecBar;
argStruct.sinusoidArgs.pitchGlideTypeVec = cell(1, Nf);
for f=1:Nf
    argStruct.sinusoidArgs.pitchGlideTypeVec{f} = 'none';
end
argStruct.sinusoidArgs.zcArgsVec = zeros(1, Nf);

% loopback FM zc parameters
argStruct.zcArgs = struct();
argStruct.zcArgs.BVec = BVec;
argStruct.zcArgs.BMat = repmat(BVec, Nf,1);
argStruct.zcArgs.BEndVec = zeros(1, Nf);
argStruct.zcArgs.gVec = zeros(1, Nf);
argStruct.zcArgs.fcVec = fcVecBar;
argStruct.zcArgs.f0EndVec = fcVecBar;
argStruct.zcArgs.BGlideTypeVec = cell(1, Nf);
for f=1:Nf
    argStruct.zcArgs.BGlideTypeVec{f} = 'useB';
end

% === marimba synthesis - traditional MS and loopback FM MS ===

% marimba - traditional MS
[marimba_s, marimba_s_Mat] = loopbackFMMS('s', env, argStruct, fs);

% marimba - loopback FM z0
[marimba_zc, marimba_zc_Mat] = loopbackFMMS('zc', env, argStruct, fs);

% === marimba synthesis - time-varying APF ===

% parameters
TVAPFParams.fpiVec = argStruct.sinusoidArgs.f0Vec;
TVAPFParams.MVec = [1000 2000 1500 3000 2300 6000 2500];
TVAPFParams.fmVec = [78, 31, 42, 83, 100, 400, 300];
TVAPFParams.fbVec = [100 200 300 400 500 600 700];

% time-varying APF - traditional MS
[yMarimba_s, yMarimba_s_Mat] = applyTimeVaryingAPF2(marimba_s_Mat, env, fs, TVAPFParams);

% time-varying APF - loopback FM z0
[yMarimba_zc, yMarimba_zc_Mat] = applyTimeVaryingAPF2(marimba_zc_Mat, env, fs, TVAPFParams);

% === marimba synthesis - commuted synthesis ===

% commuted synthesis parameters
resIRWav = '../loopbackFMPercSynth/resonatorIRs/marimbaTube.wav';
excitationType = 'rc';
excitationParams = struct();
excitationParams.winLength = 8;

% marimba - traditional MS commuted synthesis
marimba_s_CS = applyCommutedSynthesis(real(marimba_s), resIRWav, excitationType, excitationParams, fs);
%marimba_s_CS = marimba_s_CS .* wg;

% marimba - loopback FM z0 commuted synthesis
marimba_zc_CS = applyCommutedSynthesis(real(marimba_zc), resIRWav, excitationType, excitationParams, fs);
%marimba_zc_CS = marimba_zc_CS .* wg;

% marimba - traditional MS with time-varying APF - commuted synthesis
%yMarimba_s_CS = applyCommutedSynthesis(real(yMarimba_s), resIRWav, excitationType, excitationParams, fs);
%yMarimba_s_CS = yMarimba_s_CS .* wg;

% marimba - loopback FM z0 with time-varying APF - commuted synthesis
%yMarimba_zc_CS = applyCommutedSynthesis(real(yMarimba_zc), resIRWav, excitationType, excitationParams, fs);
%yMarimba_zc_CS = yMarimba_zc_CS .* wg;

% === marimba synthesis - plots ===

figure
subplot(411)
spectrogram(real(marimba_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS with CS');
%xlim([0 600])
ylim([0 10])
set(gca, 'FontSize', 15);

subplot(412)
spectrogram(real(marimba_zc_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with CS');
%xlim([0 600])
ylim([0 10])
set(gca, 'FontSize', 15);

subplot(413)
spectrogram(real(yMarimba_s), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS with time-varying APFs');
%xlim([0 600])
ylim([0 10])
set(gca, 'FontSize', 15);

subplot(414)
spectrogram(real(yMarimba_zc), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with time-varying APFs');
%xlim([0 600])
ylim([0 10])
set(gca, 'FontSize', 15);

sgtitle('Marimba synthesis example', 'FontSize', 15)

if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    saveas(gcf, [figDir filePrefix 'marimba'], 'epsc')
end

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_marimba_traditionalMS' '.wav'], ...
                scaleForSavingAudio(real(marimba_s_CS)), fs);
    audiowrite([audioDir 'hsu_marimba_loopbackFM' '.wav'], ...
                scaleForSavingAudio(real(marimba_zc_CS)), fs);
    audiowrite([audioDir 'hsu_marimba_traditionalMSTimeVaryingAPF' '.wav'], ...
                scaleForSavingAudio(real(yMarimba_s)), fs);
    audiowrite([audioDir 'hsu_marimba_loopbackFMTimeVaryingAPF' '.wav'], ...
                scaleForSavingAudio(real(yMarimba_zc)), fs);
end


%% Figure 5.11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wood block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 0.5;

% === derived parameters === 
% feedback FM
g = 0.9999;         % pitch glide coefficient
gWood = 0.9995;

N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)'; % feedback FM pitch glide coefficient

% === analyze a wood block recording to get the modal frequencies ===
[wb, fs2] = audioread('../loopbackFMPercSynth/resonatorIRs/woodBlock_forAnalysis.wav');
wb = wb(:,1);

Nwb = length(wb);
Nfft = 2^nextpow2(Nwb);
faxis = (fs2/2) * linspace(0, 1, Nfft/2+1);

WB = fft(wb, Nfft);
WBPos = WB(1:Nfft/2+1);
WBPosdB = 20*log10(abs(WBPos)/max(abs(WBPos)));

[pks, locs] = findpeaks(WBPosdB, 'minpeakheight', -60, 'minpeakdist', 80);
[spks, sinds] = sort(pks, 'descend');
spks = spks(1:11);
slocs = locs(sinds(1:11));

% plot(faxis, WBPosdB);
% hold on
% for i=1:length(spks)
%     plot(faxis(slocs(i)), spks(i), 'r*');
% end

fVecWB = faxis(slocs);
fVecWB = sort(fVecWB, 'ascend');
fVecWB = fVecWB(2:end); % drop the lowest one at like 5Hz

Nf = length(fVecWB);

% === wood block - set up decay envelopes ===

% want the lowest modal frequencies to have a larger starting amplitude 
% than the highest modal frequencies
aStart = zeros(1, Nf);
aStart(1) = 1;              % highest amplitude for lowest modal freq
aStart(Nf) = 0.001;          % lowest amplitude for highest modal freq
tau = -(Nf-1) / log(aStart(Nf));
a0 = 1/exp(-1/tau);
aStart(2:Nf-1) = a0*exp(-(2:Nf-1)/tau);     % exponentially decreasing 
                                            % starting amplitudes

e = g.^(linspace(0, N, N));   % exponential decay
eWood = gWood.^(linspace(0, N, N)); % exponential decay for woodblock

envWood = zeros(Nf, N);
for i=1:Nf
    envWood(i,:) = aStart(i) * eWood;
end

wg = envWood(1,:)';

% === wood block synthesis - loopback FM parameters ===

% traditional MS parameters
argStruct.sinusoidArgs.f0Vec = fVecWB;
argStruct.sinusoidArgs.f0EndVec = fVecWB;
argStruct.sinusoidArgs.pitchGlideTypeVec = cell(1, Nf);
for f=1:Nf
    argStruct.sinusoidArgs.pitchGlideTypeVec{f} = 'none';
end
argStruct.sinusoidArgs.zcArgsVec = zeros(1, Nf);

% loopback FM z0 parameters
argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = fVecWB;
argStruct.z0Args.f0EndVec = fVecWB / 4;
argStruct.z0Args.pitchGlideTypeVec = cell(1, Nf);
for f=1:Nf
    argStruct.z0Args.pitchGlideTypeVec{f} = 'lin';
    argStruct.z0Args.b0Mat(f,:) = linspace(0.3*((Nf-(f-1))/Nf), 0.001, N);
    argStruct.z0Args.zcArgsVec(f) = struct();
end

% === wood block synthesis - traditional MS and loopback FM MS ===

% wood block - traditional MS
[woodBlock_s, woodBlock_s_Mat] = loopbackFMMS('s', envWood, argStruct, fs);

% wood block - loopback FM z0
[woodBlock_z0, woodBlock_z0_Mat] = loopbackFMMS('z0', envWood, argStruct, fs);

% === snare synthesis - time-varying allpass filter ===

TVAPFParams.fpiVec = argStruct.sinusoidArgs.f0Vec;
TVAPFParams.fbVec = 2800 * ones(1, Nf);
TVAPFParams.MVec = 1100 * ones(1, Nf);
TVAPFParams.fmVec = 100 * ones(1, Nf);

% time-varying APF - traditional MS
[yWoodBlock_s, yWoodBlock_s_Mat] = applyTimeVaryingAPF2(woodBlock_s_Mat, envWood, fs, TVAPFParams);

% time-varying APF - loopback FM z0
[yWoodBlock_z0, yWoodBlock_z0_Mat] = applyTimeVaryingAPF2(woodBlock_z0_Mat, envWood, fs, TVAPFParams);

% === woodBlock synthesis - commuted synthesis ===

% commuted synthesis parameters
resIRWav = '../loopbackFMPercSynth/resonatorIRs/113620__vidsyn__miscsoftnaturalgtrloud2-2.wav';
excitationType = 'nb';
excitationParams = struct();
excitationParams.durNB = 0.2;%0.001;
excitationParams.lowFreq = 120;
excitationParams.highFreq = 12000;

% woodBlock - traditional MS commuted synthesis
woodBlock_s_CS = applyCommutedSynthesis(real(woodBlock_s), resIRWav, excitationType, excitationParams, fs);
woodBlock_s_CS = woodBlock_s_CS .* wg;

if size(real(woodBlock_s)') == size(wg)
    woodBlock_s = real(woodBlock_s)' .* wg;
end

% woodBlock - loopback FM z0 commuted synthesis
woodBlock_z0_CS = applyCommutedSynthesis(real(woodBlock_z0), resIRWav, excitationType, excitationParams, fs);
woodBlock_z0_CS = woodBlock_z0_CS .* wg;

if size(real(woodBlock_z0)') == size(wg)
    woodBlock_z0 = real(woodBlock_z0)' .* wg;
end

% woodBlock - traditional MS with time-varying APF - commuted synthesis
%yWoodBlock_s_CS = applyCommutedSynthesis(real(yWoodBlock_s), resIRWav, excitationType, excitationParams, fs);
%yWoodBlock_s = yWoodBlock_s' .* wg;

% woodBlock - loopback FM z0 with time-varying APF - commuted synthesis
%yWoodBlock_z0_CS = applyCommutedSynthesis(real(yWoodBlock_z0), resIRWav, excitationType, excitationParams, fs);
%yWoodBlock_z0 = yWoodBlock_z0' .* wg;

% === woodBlock synthesis - plot ===
figure
subplot(411)
spectrogram(real(woodBlock_s), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS');
xlim([0 500])
ylim([0 9])
set(gca, 'FontSize', 15);

subplot(412)
spectrogram(real(woodBlock_z0), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS');
xlim([0 500])
ylim([0 9])
set(gca, 'FontSize', 15);

subplot(413)
spectrogram(real(woodBlock_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS with CS');
xlim([0 500])
ylim([0 9])
set(gca, 'FontSize', 15);

subplot(414)
spectrogram(real(woodBlock_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with CS');
xlim([0 500])
ylim([0 9])
set(gca, 'FontSize', 15);

sgtitle('Wood block synthesis example', 'FontSize', 15)

if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    saveas(gcf, [figDir filePrefix 'woodBlock'], 'epsc')
end

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_woodBlock_traditionalMS' '.wav'], ...
                scaleForSavingAudio(real(woodBlock_s)), fs);
    audiowrite([audioDir 'hsu_woodBlock_loopbackFM' '.wav'], ...
                scaleForSavingAudio(real(woodBlock_z0)), fs);
    audiowrite([audioDir 'hsu_woodBlock_traditionalMSCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(woodBlock_s_CS)), fs);
    audiowrite([audioDir 'hsu_woodBlock_loopbackFMCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(woodBlock_z0_CS)), fs);
end

%% Figure 5.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tom tom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 2.0;

% === tom tom - derived parameters ===

% feedback FM
B = 0.9;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient

% === tom tom - ideal membrane modal frequencies ===
% from Science of Percussion Instruments

% TOM-TOM (page 26)
f1 = 142;
fVecMembrane = f1 * [1, 2.15, 3.17, 3.42, 4.09, 4.80, 4.94];

Nf = length(fVecMembrane);

% center frequencies for feedback FM if sounding frequencies = membrane
% modal freqeuncies
fcVecMembrane = fVecMembrane./(sqrt(1-B^2));

% === tom tom - set up decay envelopes ===
% want the lowest modal frequencies to have a larger starting amplitude 
% than the highest modal frequencies
aStart = zeros(1, Nf);
aStart(1) = 1;              % highest amplitude for lowest modal freq
aStart(Nf) = 0.01;          % lowest amplitude for highest modal freq
tau = -(Nf-1) / log(aStart(Nf));
a0 = 1/exp(-1/tau);
aStart(2:Nf-1) = a0*exp(-(2:Nf-1)/tau);     % exponentially decreasing 
                                            % starting amplitudes

e = g.^(linspace(0, N, N));   % exponential decay

env = zeros(Nf, N);
for i=1:Nf
    env(i,:) = aStart(i) * e;
end

wg = env(1,:)';

% === tom tom synthesis - loopback FM parameters ===
argStruct = struct();

% traditional MS parameters
argStruct.sinusoidArgs.f0Vec = fVecMembrane;
argStruct.sinusoidArgs.f0EndVec = fVecMembrane;
argStruct.sinusoidArgs.pitchGlideTypeVec = cell(1, Nf);
for f=1:Nf
    argStruct.sinusoidArgs.pitchGlideTypeVec{f} = 'none';
end
argStruct.sinusoidArgs.zcArgsVec = zeros(1, Nf);

% loopback FM z0 parameters
argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = fVecMembrane;
argStruct.z0Args.f0EndVec = fVecMembrane / 1.3;
argStruct.z0Args.pitchGlideTypeVec = cell(1, Nf);
b0 = -0.3;

for f=1:Nf
    argStruct.z0Args.pitchGlideTypeVec{f} = 'exp';
    argStruct.z0Args.b0Mat(f,:) = b0*ones(1, N);
    argStruct.z0Args.zcArgsVec(f) = struct();
end
for f=1:Nf % not sure why i have to do this here...
    argStruct.z0Args.zcArgsVec(f).T60 = 0.1;
end


% === tom tom synthesis - traditional MS and loopback FM MS ===

% tom tom - traditional MS
[tomtom_s, tomtom_s_Mat] = loopbackFMMS('s', env, argStruct, fs);

% tom tom - loopback FM z0
[tomtom_z0, tomtom_z0_Mat] = loopbackFMMS('z0', env, argStruct, fs);


% === tom tom synthesis - commuted synthesis ===

% commuted synthesis parameters
resIRWav = '../loopbackFMPercSynth/resonatorIRs/taiko/taiko2_fixedOffset.wav';
excitationType = 'nb';
excitationParams = struct();
excitationParams.durNB = 0.02;
excitationParams.lowFreq = 120;
excitationParams.highFreq = 4000;

% tom tom - traditional MS commuted synthesis
tomtom_s_CS = applyCommutedSynthesis(real(tomtom_s), resIRWav, excitationType, excitationParams, fs);
tomtom_s_CS = tomtom_s_CS(1:N);
%tomtom_s_CS = tomtom_s_CS .* wg;

% tom tom - loopback FM z0 commuted synthesis
tomtom_z0_CS = applyCommutedSynthesis(real(tomtom_z0), resIRWav, excitationType, excitationParams, fs);
tomtom_z0_CS = tomtom_z0_CS(1:N);
%tomtom_z0_CS = tomtom_z0_CS .* wg;


% === tom tom synthesis - plot ===
figure
subplot(411)
spectrogram(real(tomtom_s), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS');
xlim([0 dur])
ylim([0 7])
set(gca, 'FontSize', 15);

subplot(412)
spectrogram(real(tomtom_z0), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS');
xlim([0 dur])
ylim([0 7])
set(gca, 'FontSize', 15);

subplot(413)
spectrogram(real(tomtom_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS with CS');
xlim([0 dur])
ylim([0 7])
set(gca, 'FontSize', 15);

subplot(414)
spectrogram(real(tomtom_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with CS');
xlim([0 dur])
ylim([0 7])
set(gca, 'FontSize', 15);

sgtitle('Tom tom drum synthesis example', 'FontSize', 15)

if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    saveas(gcf, [figDir filePrefix 'tomtom'], 'epsc')
end

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_tomtom_traditionalMS' '.wav'], ...
                scaleForSavingAudio(real(tomtom_s)), fs);
    audiowrite([audioDir 'hsu_tomtom_loopbackFM' '.wav'], ...
                scaleForSavingAudio(real(tomtom_z0)), fs);
    audiowrite([audioDir 'hsu_tomtom_traditionalMSCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(tomtom_s_CS)), fs);
    audiowrite([audioDir 'hsu_tomtom_loopbackFMCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(tomtom_z0_CS)), fs);
end


%% Figure 5.13, 5.14, and 5.15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circular Plate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 2;

% === circular plate - derived parameters ===

% feedback FM
B = 0.99;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient

% === circular plate - modal frequencies with a clamped edge ===
% from Science of Percussion Instruments (page 80)

h = 0.005;
a = 0.09;

E = 2*10^11;
v = 0.3;
rho = 7860;
cL = sqrt(E/(rho*(1 - v^2)));

%f01 = 0.4694*cL*h/a^2;
%f01 = 2000;
%fVecCP = f01 * [1, 2.08, 3.41, 5, 6.82, 3.89, 5.95, 8.28, 10.87, 13.71, 8.72, 11.75, 15.06, 18.63, 22.47];
% SIMPLY-SUPPORTED CIRCULAR PLATE

f01 = 0.2287 * cL * h / a^2;
fVecCP = f01 * [1, 2.80, 5.15, 5.98, 9.75, 14.09, 14.91, 20.66, 26.99];

Nf = length(fVecCP);

% center frequencies for feedback FM if sounding frequencies = membrane
% modal freqeuncies
fcVecCP = fVecCP./(sqrt(1-B^2));

% === circular plate - set up decay envelopes ===
% want the lowest modal frequencies to have a larger starting amplitude 
% than the highest modal frequencies
aStart = zeros(1, Nf);
aStart(1) = 1;              % highest amplitude for lowest modal freq
aStart(Nf) = 0.5;          % lowest amplitude for highest modal freq
tau = -(Nf-1) / log(aStart(Nf));
a0 = 1/exp(-1/tau);
aStart(2:Nf-1) = a0*exp(-(2:Nf-1)/tau);     % exponentially decreasing 
                                            % starting amplitudes

e = g.^(linspace(0, N, N));   % exponential decay

env = zeros(Nf, N);
for i=1:Nf
    env(i,:) = aStart(i) * e;
end

% === circular plate synthesis - loopback FM parameters ===

% traditional MS parameters
argStruct.sinusoidArgs.f0Vec = fVecCP;
argStruct.sinusoidArgs.f0EndVec = fVecCP;
argStruct.sinusoidArgs.pitchGlideTypeVec = cell(1, Nf);
for f=1:Nf
    argStruct.sinusoidArgs.pitchGlideTypeVec{f} = 'none';
end
argStruct.sinusoidArgs.zcArgsVec = zeros(1, Nf);

% loopback FM zc parameters
argStruct.zcArgs = struct();
BVec = linspace(0.91, 0.9, N);
argStruct.zcArgs.BVec = BVec;
argStruct.zcArgs.BMat = repmat(BVec, Nf,1);
argStruct.zcArgs.BEndVec = zeros(1, Nf);
argStruct.zcArgs.gVec = zeros(1, Nf);
argStruct.zcArgs.fcVec = fcVecCP;
argStruct.zcArgs.f0EndVec = fcVecCP;
argStruct.zcArgs.BGlideTypeVec = cell(1, Nf);
for f=1:Nf
    argStruct.zcArgs.BGlideTypeVec{f} = 'useB';
end

% === circular plate synthesis - traditional MS and loopback FM MS ===

% marimba - traditional MS
[circularPlate_s, circularPlate_s_Mat] = loopbackFMMS('s', env, argStruct, fs);

% marimba - loopback FM z0
[circularPlate_zc, circularPlate_zc_Mat] = loopbackFMMS('zc', env, argStruct, fs);

% === circular plate synthesis - time-varying APF ===

% parameters
TVAPFParams.fpiVec = argStruct.sinusoidArgs.f0Vec;
TVAPFParams.MVec = 11500*ones(1, Nf);
TVAPFParams.fmVec = 1000*ones(1, Nf);
TVAPFParams.fbVec = 2750*ones(1, Nf);

% time-varying APF - traditional MS
[yCircularPlate_s, yCircularPlate_s_Mat] = applyTimeVaryingAPF2(circularPlate_s_Mat, env, fs, TVAPFParams);

% time-varying APF - loopback FM z0
[yCircularPlate_zc, yCircularPlate_zc_Mat] = applyTimeVaryingAPF2(circularPlate_zc_Mat, env, fs, TVAPFParams);

% === circular plate synthesis - commuted synthesis ===

% commuted synthesis parameters
resIRWav = '../loopbackFMPercSynth/resonatorIRs/CarpenterCenter.wav';
excitationType = 'rc';
excitationParams = struct();
excitationParams.winLength = 16;

% circular plate - traditional MS commuted synthesis
circularPlate_s_CS = applyCommutedSynthesis(real(circularPlate_s), resIRWav, excitationType, excitationParams, fs);

% circular plate - loopback FM z0 commuted synthesis
circularPlate_zc_CS = applyCommutedSynthesis(real(circularPlate_zc), resIRWav, excitationType, excitationParams, fs);

% circular plate - traditional MS commuted synthesis
yCircularPlate_s_CS = applyCommutedSynthesis(real(yCircularPlate_s), resIRWav, excitationType, excitationParams, fs);

% circular plate - loopback FM z0 commuted synthesis
yCircularPlate_zc_CS = applyCommutedSynthesis(real(yCircularPlate_zc), resIRWav, excitationType, excitationParams, fs);

% === circular plate synthesis - plots ===

figure
subplot(211)
spectrogram(real(circularPlate_s), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS');
%xlim([0 600])
ylim([0 15])
set(gca, 'FontSize', 15);

subplot(212)
spectrogram(real(circularPlate_zc), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS');
%xlim([0 600])
ylim([0 15])
set(gca, 'FontSize', 15);

sgtitle('Circular plate synthesis example', 'FontSize', 15)
if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 5];
    saveas(gcf, [figDir filePrefix 'circularPlate'], 'epsc')
end

figure
subplot(211)
spectrogram(real(circularPlate_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS with CS');
%xlim([0 600])
ylim([0 15])
set(gca, 'FontSize', 15);

subplot(212)
spectrogram(real(circularPlate_zc_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with CS');
%xlim([0 600])
ylim([0 15])
set(gca, 'FontSize', 15);

sgtitle('Circular plate synthesis example with CS', 'FontSize', 15)

if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 5];
    saveas(gcf, [figDir filePrefix 'circularPlateCS'], 'epsc')
end

figure
subplot(211)
spectrogram(real(yCircularPlate_s), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS');
%xlim([0 600])
ylim([0 15])
set(gca, 'FontSize', 15);

subplot(212)
spectrogram(real(yCircularPlate_zc), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS');
%xlim([0 600])
ylim([0 15])
set(gca, 'FontSize', 15);

sgtitle('Circular plate synthesis w/ time-varying APFs and CS', 'FontSize', 15)

if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 5];
    saveas(gcf, [figDir filePrefix 'circularPlateTimeVaryingAPFsCS'], 'epsc')
end

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_circularPlate_traditionalMS' '.wav'], ...
                scaleForSavingAudio(real(circularPlate_s)), fs);
    audiowrite([audioDir 'hsu_circularPlate_loopbackFM' '.wav'], ...
                scaleForSavingAudio(real(circularPlate_zc)), fs);
    audiowrite([audioDir 'hsu_circularPlate_traditionalMSCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(circularPlate_s_CS)), fs);
    audiowrite([audioDir 'hsu_circularPlate_loopbackFMCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(circularPlate_zc_CS)), fs);
    audiowrite([audioDir 'hsu_circularPlate_traditionalMSTimeVaryingAPFCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(yCircularPlate_s_CS)), fs);
    audiowrite([audioDir 'hsu_circularPlate_loopbackFMTimeVaryingAPFCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(yCircularPlate_zc_CS)), fs);
end