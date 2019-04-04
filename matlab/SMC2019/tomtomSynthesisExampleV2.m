% tomtomSynthesisExampleV2.m
%
% Using the tom tom modal frequencies, this script performs:
%   a. modal synthesis
%   b. the time-varying allpass filtering for sinusoids used in 
%       traditional modal synthesis
%   c. loopback FM modal synthesis
%   d. time-varying allpass filtering on loopback FM oscillators
% The tom tom modal synthesis signals are saved to
%   timeVaryingAPF/audioExamples/tomtom directory. The file
%   tomtomCommutedSynthesisAndPlot.m can be called to perform 
%   commuted synthesis on these saved audio files and to plot the file
%   results.

addpath(genpath('../proofOfConcept'));

%% input parameters

% general
fs = 44100;
dur = 1;
plotSpectrograms = 1;
writeAudioFiles = 1;
outDir = 'audioExamples/tomtom/V2/';

if ~exist(outDir)
    mkdir(outDir)
end

%% derived parameters

% feedback FM
B = 0.9;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient


%% ideal membrane modal frequencies
% from Science of Percussion Instruments

% TOM-TOM (page 26)
f1 = 142;
fVecMembrane = f1 * [1, 2.15, 3.17, 3.42, 4.09, 4.80, 4.94];

Nf = length(fVecMembrane);

% center frequencies for feedback FM if sounding frequencies = membrane
% modal freqeuncies
fcVecMembrane = fVecMembrane./(sqrt(1-B^2));


%% set up decay envelopes
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

%% modal synthesis (no pitch glide)
[yMS, yMSMat] = modalSynthSine(fVecMembrane', env, fs, 0);


%% loopback FM using stretched APF formulation

%fVecMembraneStart = fVecMembrane*2;
%[yLBFM, ~] = stretchedAPFSynthesis(fVecMembraneStart, -0.3, env, fs, fVecMembrane, 'expB');

%fVecMembraneEnd = fVecMembrane/2;
%[yLBFM, ~] = stretchedAPFSynthesis(fVecMembrane, -0.3, env, fs, fVecMembraneEnd, 'exp');

fVecMembraneEnd = fVecMembrane/2;
b0 = -0.3;
[yLBFM, ~] = stretchedAPFSynthesis(fVecMembrane, b0, env, fs, fVecMembraneEnd, 'linearB');


%% find the pitch glide that I am using in my examples
%f0VecBar2 = fcVecBar2(:) .* sqrt(1 - BVec2.^2);

% 'exp' type pitch glide with starting frequency at fVecMembrane_w0
% g = 0.9999;
% BB = g.^(0:N-1);
% f0Vec = fVecMembrane' .* sqrt(1 - BB.^2);

% nT = 0:T:(dur-T);
% tau = 1./(log(fVecMembraneEnd./fVecMembrane));
% f0Vec = fVecMembrane(:) .* exp(nT./tau(:));
% Theta = 2*pi*fVecMembrane(:) .* tau(:) .* exp(nT./tau(:));


B1 = -2*b0/(b0^2 + 1);
fc = fVecMembrane(1)/(sqrt(1 - B1^2));
B2 = sqrt(1 - (fVecMembraneEnd(1)/fc)^2);

k = (B2 - B1)/(N - 1);
l = B1 - k;
n = 0:N-1;
BVec2 = k*n + l;

fcVecMembrane2(:) = fVecMembrane(:)/(sqrt(1 - B1^2));
f0Vec = fcVecMembrane2(:) * sqrt(1 - BVec2.^2);

Theta = (fcVecMembrane2(:)*T/(2*k)) .* ((k*n + l) .* sqrt(1 - (k*n + l).^2) + asin(k*n + l));


% check with a plot
if plotSpectrograms == 1
    figure
    subplot(211)
    spectrogram(real(yLBFM), hann(256), 128, 1024, fs, 'yaxis');
    title('loopback FM MS synthesis spectrogram')
    hold on
    plot(linspace(0, 500, N), f0Vec/1000, 'r')
end

%% modal synthesis (sinusoidal) with exponential decay using the modes
% u = sqrt(1 - BB.^2);
% C = 0;                      % constant of integration
% Theta = (2*pi*fVecMembrane(:)*T/log(g)) .* (u - atanh(u) + C);
% f0VecPG = Theta/(2*pi);
%         
% [yMSPG, yMSMat] = modalSynthSine(f0VecPG, env, fs, 1);
[yMSPG, yMSMat] = modalSynthSine(Theta, env, fs, 1);

%% modal synthesis (bandpass) with exponential decay using the modes
[yMSBP, yMSBPMat] = modalSynthBP(f0Vec, env, fs);

if plotSpectrograms == 1
    figure
    subplot(311)
    spectrogram(real(yMSPG), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, 1000, N), f0Vec/1000, 'r');
    ylim([0 3])
    title('MS spectrogram (sinusoid)')
    subplot(312)
    spectrogram(real(yMSBP), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, 1000, N), f0Vec/1000, 'r');
    ylim([0 3])
    title('MS spectrogram (bandpass)')
    subplot(313)
    spectrogram(real(yLBFM), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, 1000, N), f0Vec/1000, 'r');
    ylim([0 3])
    title('loopback FM MS spectrogram')
end


%% Other loopback FM synthesis for the sound example webpage

BVec2 = g.^(0:N-1)';
[yLBFMMemb3, yLBFMMembMat3] = feedbackFMSynthesis(fcVecMembrane, BVec2, env, fs);

% find the pitch glide that this is using
f0VecMembrane = fcVecMembrane .* sqrt(1 - BVec2.^2);

fVecEnd = f0VecMembrane(end,:);
% modal synthesis (with a pitch glide)
u = sqrt(1 - BVec2.^2)';
C = 0;                      % constant of integration
Theta = (2*pi*fVecEnd(:).*T/log(g)) .* (u - atanh(u) + C);
[yLBFMMemb3MSPG, ~] = modalSynthSine(Theta/(2*pi), env, fs, 1);


%% time-varying allpass traditional modal synthesis
% [for dissertation, not SMC]
% 
% % parameters
% TVAPFParams.M = 11500;
% TVAPFParams.f_m = 100;
% TVAPFParams.f_b = 2750;
% 
% [yTVAPFMS, yTVAPFWBMat2, TVAPFParams1] = TVAPFSynthesis(fVecMembrane, env, TVAPFParams, 0, fs);
% 
% if plotSpectrograms == 1
%     figure
%     subplot(211)
%     spectrogram(real(yMS), hann(256), 128, 1024, fs, 'yaxis');
%     title('modal synthesis spectrogram')
%     subplot(212)
%     spectrogram(real(yTVAPFMS), hann(256), 128, 1024, fs, 'yaxis');
%     title('time-varying APF modal synthesis spectrogram')
% end


%% time-varying allpass loopback FM (stretched APF)
% [for dissertation, not SMC]

% parameters
% TVAPFParams.M = 1500;
% TVAPFParams.f_m = 100;
% TVAPFParams.f_b = 2750;
% 
% % loopback FM as stretched APF
% [yLBFM, ~] = stretchedAPFSynthesis(fVecMembrane_w0, -0.9, env, fs, fVecMembrane_w0/2, 'exp');
% 
% % time-varying allpass filter each loopback FM oscillator
% yTVAPFLB = zeros(1, N);
% envOnes = ones(size(env));
% [~, yLBFMMat] = stretchedAPFSynthesis(fVecMembrane_w0, -0.9, envOnes, fs, fVecMembrane_w0/2, 'exp');
% for i=1:Nf
%     f_pi = fVecMembrane(i);
%     %yTVAPFLB = yTVAPFLB + TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams1.f_b(i), TVAPFParams1.M(i), TVAPFParams1.f_m(i), fs);
%     y1 = TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams.f_b, TVAPFParams.M, TVAPFParams.f_m, fs);
%     yTVAPFLB = yTVAPFLB + y1 .* env(i,:);
% end
% 
% if plotSpectrograms == 1
%     figure
%     subplot(211)
%     spectrogram(real(yLBFM), hann(256), 128, 1024, fs, 'yaxis');
%     title('loopback FM MS spectrogram')
%     subplot(212)
%     spectrogram(real(yTVAPFLB), hann(256), 128, 1024, fs, 'yaxis');
%     title('time-varying APF loopback FM MS spectrogram')
% end


%% Write audio files

if writeAudioFiles == 1
    audiowrite([outDir 'modalSynthesis.wav'], scaleForSavingAudio(real(yMS)), fs);
    audiowrite([outDir 'modalSynthesisPitchGlide.wav'], scaleForSavingAudio(real(yMSPG)), fs);
    %audiowrite([outDir 'TVAPFModalSynthesis.wav'], scaleForSavingAudio(real(yTVAPFMS)), fs);
    audiowrite([outDir 'yLBFMMemb3ModalSynthesisPitchGlide.wav'], scaleForSavingAudio(real(yLBFMMemb3MSPG)), fs);
    
    audiowrite([outDir 'loopbackFM.wav'], scaleForSavingAudio(real(yLBFM)), fs);
    audiowrite([outDir 'yLBFMMemb3.wav'], scaleForSavingAudio(real(yLBFMMemb3)), fs);
    %audiowrite([outDir 'TVAPFLoopbackFM.wav'], scaleForSavingAudio(real(yTVAPFLB)), fs);
end
