% reverse reverb circular plate
% hihat

addpath(genpath('../loopbackFMPercSynth/'));
addpath(genpath('../helperFunctions/'));
savePlots = 1;
saveAudio = 1;

figDir = 'figures/';
audioDir = 'audioExamples/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circular Plate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 3.0;

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

%e = g.^(linspace(0, N, N));   % exponential decay

T60 = 2;
tau  = -T60/log(0.001);
e = 1 * exp(-(0:N-1)*T/tau);
e = fliplr(e);

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
BVec = linspace(0.99, 0.99, N);
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

% circular plate - traditional MS
[circularPlate_s, circularPlate_s_Mat] = loopbackFMMS('s', env, argStruct, fs);

% circular plate - loopback FM zc
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
resIRWav = '../loopbackFMPercSynth/resonatorIRs/reverseCarpenterCenter.wav';
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

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_reverseCP_traditionalMS' '.wav'], ...
                scaleForSavingAudio(real(circularPlate_s)), fs);
    audiowrite([audioDir 'hsu_reverseCP_loopbackFM' '.wav'], ...
                scaleForSavingAudio(real(circularPlate_zc)), fs);
    audiowrite([audioDir 'hsu_reverseCP_traditionalMSCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(circularPlate_s_CS)), fs);
    audiowrite([audioDir 'hsu_reverseCP_loopbackFMCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(circularPlate_zc_CS)), fs);
    audiowrite([audioDir 'hsu_reverseCP_traditionalMSTimeVaryingAPFCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(yCircularPlate_s_CS)), fs);
    audiowrite([audioDir 'hsu_reverseCP_loopbackFMTimeVaryingAPFCommutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(yCircularPlate_zc_CS)), fs);
end