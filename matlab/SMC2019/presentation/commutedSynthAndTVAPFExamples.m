% commuted synthesis and time-varying APFs examples

addpath(genpath('../../loopbackFMPercSynth/'));
addpath(genpath('../../helperFunctions/'));

% directories for saving
figDir = 'figures/';
audioDir = 'audioExamples/';

savePlots = 1;
saveAudio = 1;

%% MARIMBA - COMMUTED SYNTHESIS EXAMPLE
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

% marimba - loopback FM z0
[marimba_zc, marimba_zc_Mat] = loopbackFMMS('zc', env, argStruct, fs);

% === marimba synthesis - commuted synthesis ===

% commuted synthesis parameters
resIRWav = '../../loopbackFMPercSynth/resonatorIRs/CarpenterCenter.wav';
excitationType = 'rc';
excitationParams = struct();
excitationParams.winLength = 8;

% marimba - loopback FM z0 commuted synthesis
marimba_zc_CS = applyCommutedSynthesis(real(marimba_zc), resIRWav, excitationType, excitationParams, fs);


% === marimba synthesis - plots ===

figure

subplot(211)
spectrogram(real(marimba_zc), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS');
%xlim([0 600])
ylim([0 12])
set(gca, 'FontSize', 15);

subplot(212)
spectrogram(real(marimba_zc_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with CS');
%xlim([0 600])
ylim([0 12])
set(gca, 'FontSize', 15);

sgtitle('Marimba synthesis example', 'FontSize', 15)

if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 5];
    saveas(gcf, [figDir 'marimba_commutedSynthesis'], 'epsc')
end

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_marimba_loopbackFM' '.wav'], ...
                scaleForSavingAudio(real(marimba_zc)), fs);
    audiowrite([audioDir 'hsu_marimba_loopbackFM_commutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(marimba_zc_CS)), fs);
end


%% SNARE DRUM - TIME-VARYING APF EXAMPLE
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

% snare - loopback FM z0
[snare_z0, snare_z0_Mat] = loopbackFMMS('z0', env, argStruct, fs);

% === snare synthesis - time-varying allpass filter ===

TVAPFParams.fpiVec = argStruct.sinusoidArgs.f0Vec;
TVAPFParams.fbVec = [100 200];
TVAPFParams.MVec = [1000 250];
TVAPFParams.fmVec = [300 2000];


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

% snare - loopback FM z0 commuted synthesis
snare_z0_CS = applyCommutedSynthesis(real(snare_z0), resIRWav, excitationType, excitationParams, fs);
snare_z0_CS = snare_z0_CS .* wg;

% snare - loopback FM z0 with time-varying APF - commuted synthesis
ySnare_z0_CS = applyCommutedSynthesis(real(ySnare_z0), resIRWav, excitationType, excitationParams, fs);
ySnare_z0_CS = ySnare_z0_CS .* wg;

% === snare synthesis - plot ===
figure

subplot(211)
spectrogram(real(snare_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with CS');
xlim([0 600])
ylim([0 15])
set(gca, 'FontSize', 15);

subplot(212)
spectrogram(real(ySnare_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS w/ CS and time-varying APFs');
xlim([0 600])
ylim([0 15])
set(gca, 'FontSize', 15);

sgtitle('Snare drum synthesis example', 'FontSize', 15)

if savePlots
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 5];
    saveas(gcf, [figDir 'snareAPF'], 'epsc')
end

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_snare_loopbackFM' '.wav'], ...
                scaleForSavingAudio(real(snare_z0_CS)), fs);
    audiowrite([audioDir 'hsu_snare_loopbackFMTimeVaryingAPF' '.wav'], ...
                scaleForSavingAudio(real(ySnare_z0_CS)), fs);
end
