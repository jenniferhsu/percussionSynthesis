% tomtomExamples.m

addpath(genpath('../../loopbackFMPercSynth/'));
addpath(genpath('../../helperFunctions/'));
savePlots = 1;
saveAudio = 1;

figDir = 'figures/';
audioDir = 'audioExamples/';


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
resIRWav = '../../loopbackFMPercSynth/resonatorIRs/taiko/taiko2_fixedOffset.wav';
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
subplot(221)
spectrogram(real(tomtom_s), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS');
xlim([0 dur])
ylim([0 7])
set(gca, 'FontSize', 15);

subplot(222)
spectrogram(real(tomtom_z0), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS');
xlim([0 dur])
ylim([0 7])
set(gca, 'FontSize', 15);

subplot(223)
spectrogram(real(tomtom_s_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Traditional MS with CS');
xlim([0 dur])
ylim([0 7])
set(gca, 'FontSize', 15);

subplot(224)
spectrogram(real(tomtom_z0_CS), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Loopback FM MS with CS');
xlim([0 dur])
ylim([0 7])
set(gca, 'FontSize', 15);

sgtitle('Tom tom', 'FontSize', 15)

if savePlots
    saveas(gcf, [figDir 'tomtomExample'], 'epsc')
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
