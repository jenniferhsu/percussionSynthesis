% snareExamples.m

addpath(genpath('../../loopbackFMPercSynth/'));
addpath(genpath('../../helperFunctions/'));
savePlots = 1;
saveAudio = 1;

figDir = 'figures/';
audioDir = 'audioExamples/';

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
resIRWav = '../../loopbackFMPercSynth/resonatorIRs/impulse.wav';
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
subplot(221)
spectrogram(real(snare_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('MS with CS');
xlim([0 600])
ylim([0 5])
set(gca, 'FontSize', 10);

subplot(222)
spectrogram(real(snare_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with CS');
xlim([0 600])
ylim([0 5])
set(gca, 'FontSize', 10);

subplot(223)
spectrogram(real(ySnare_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('MS with CS and time-varying APFs');
xlim([0 600])
ylim([0 5])
set(gca, 'FontSize', 10);

subplot(224)
spectrogram(real(ySnare_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS w/ CS and time-varying APFs');
xlim([0 600])
ylim([0 5])
set(gca, 'FontSize', 10);

sgtitle('Snare drum', 'FontSize', 15)

if savePlots
    saveas(gcf, [figDir 'snareExample'], 'epsc')
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