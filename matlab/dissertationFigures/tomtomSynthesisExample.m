% tomtomSynthesisExample.m
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
dur = 0.5;
plotSpectrograms = 1;
writeAudioFiles = 1;
outDir = 'audioExamples/tomtom/';

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


%% modal synthesis with exponential decay using the modes

yMS = zeros(1, N);
nVec = 0:T:(dur-T);

for i=1:Nf
    f = fVecMembrane(i);
    yMS = yMS + (exp(1j*2*pi*f*nVec) .* env(i,:));
end

%% time-varying allpass traditional modal synthesis

% parameters
TVAPFParams.M = 100;
TVAPFParams.f_m = 250;
TVAPFParams.f_b = 900;

% TVAPFParams.M = fs/40;
% TVAPFParams.f_m = 100;
% TVAPFParams.f_b = fs/16;

[yTVAPFMS, yTVAPFWBMat2, TVAPFParams1] = TVAPFSynthesis(fVecMembrane, env, TVAPFParams, 0, fs);

if plotSpectrograms == 1
    figure
    subplot(211)
    spectrogram(real(yMS), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
    subplot(212)
    spectrogram(real(yTVAPFMS), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying APF modal synthesis spectrogram')
end


%% time-varying allpass loopback FM (stretched APF)

% parameters
TVAPFParams.M = 100;
TVAPFParams.f_m = 250;
TVAPFParams.f_b = 900;

b0 = (sqrt(1-B^2) - 1)/B;
fVecMembrane_w0 = fVecMembrane*sqrt(1 - B^2);

% loopback FM as stretched APF
[yLBFM, ~] = stretchedAPFSynthesis(fVecMembrane_w0, -0.9, env, fs, fVecMembrane_w0/2, 'exp');

% time-varying allpass filter each loopback FM oscillator
yTVAPFLB = zeros(1, N);
envOnes = ones(size(env));
[~, yLBFMMat] = stretchedAPFSynthesis(fVecMembrane_w0, -0.9, envOnes, fs, fVecMembrane_w0/2, 'exp');
for i=1:Nf
    f_pi = fVecMembrane(i);
    %yTVAPFLB = yTVAPFLB + TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams1.f_b(i), TVAPFParams1.M(i), TVAPFParams1.f_m(i), fs);
    y1 = TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams.f_b, TVAPFParams.M, TVAPFParams.f_m, fs);
    yTVAPFLB = yTVAPFLB + y1 .* env(i,:);
end

if plotSpectrograms == 1
    figure
    subplot(211)
    spectrogram(real(yLBFM), hann(256), 128, 1024, fs, 'yaxis');
    title('loopback FM MS spectrogram')
    subplot(212)
    spectrogram(real(yTVAPFLB), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying APF loopback FM MS spectrogram')
end

%% Write audio files

if writeAudioFiles == 1
    audiowrite([outDir 'modalSynthesis.wav'], scaleForSavingAudio(real(yMS)), fs);
    audiowrite([outDir 'TVAPFModalSynthesis.wav'], scaleForSavingAudio(real(yTVAPFMS)), fs);
    
    audiowrite([outDir 'loopbackFM.wav'], scaleForSavingAudio(real(yLBFM)), fs);
    audiowrite([outDir 'TVAPFLoopbackFM.wav'], scaleForSavingAudio(real(yTVAPFLB)), fs);
end
