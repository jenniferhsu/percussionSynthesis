addpath(genpath('../loopbackFMPercSynth/'));
addpath(genpath('../helperFunctions/'));
savePlots = 0;
saveAudio = 1;

figDir = 'figures/';
audioDir = 'audioExamples/';

% WOOD BLOCK TUNED TO C

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

% tune to key of C
fVecWB = fVecWB/fVecWB(1);
%fVecWB = fVecWB * 261.63;

% tune to the key of G
fVecWB = fVecWB * 196;

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


% save audio
if saveAudio
%     audiowrite([audioDir 'hsu_woodBlock_traditionalMS' '.wav'], ...
%                 scaleForSavingAudio(real(woodBlock_s)), fs);
%     audiowrite([audioDir 'hsu_woodBlock_loopbackFM' '.wav'], ...
%                 scaleForSavingAudio(real(woodBlock_z0)), fs);
%     audiowrite([audioDir 'hsu_woodBlock_traditionalMSCommutedSynthesis' '.wav'], ...
%                 scaleForSavingAudio(real(woodBlock_s_CS)), fs);

%    audiowrite([audioDir 'hsu_woodBlock_loopbackFMCommutedSynthesis_KeyC' '.wav'], ...
%                scaleForSavingAudio(real(woodBlock_z0_CS)), fs);
            
    audiowrite([audioDir 'hsu_woodBlock_loopbackFMCommutedSynthesis_KeyG' '.wav'], ...
                scaleForSavingAudio(real(woodBlock_z0_CS)), fs);
end
