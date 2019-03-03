% circularPlateSynthesisExample.m
%
% Using the circular plate modal frequencies, this script performs:
%   a. modal synthesis
%   b. the time-varying allpass filtering for sinusoids used in 
%       traditional modal synthesis
%   c. loopback FM modal synthesis
%   d. time-varying allpass filtering on loopback FM oscillators
% The circular plate  modal synthesis signals are saved to
%   timeVaryingAPF/audioExamples/circularPlate/ directory. The file
%   circularPlateCommutedSynthesisAndPlot.m can be called to perform 
%   commuted synthesis on these saved audio files and to plot the file
%   results.

addpath(genpath('../proofOfConcept'));

%% input parameters

% general
fs = 44100;
dur = 0.5;
plotSpectrograms = 1;
writeAudioFiles = 1;
outDir = 'audioExamples/circularPlate/';

if ~exist(outDir)
    mkdir(outDir)
end

%% derived parameters

% feedback FM
B = 0.99;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient


%% circular plate with a clamped edge
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


%% set up decay envelopes
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


%% modal synthesis with exponential decay using the modes

yMS = zeros(1, N);
nVec = 0:T:(dur-T);

for i=1:Nf
    f = fVecCP(i);
    yMS = yMS + (exp(1j*2*pi*f*nVec) .* env(i,:));
end


%% time-varying allpass traditional modal synthesis

% parameters
TVAPFParams.M = 11500;
TVAPFParams.f_m = 1000;
TVAPFParams.f_b = 2750;

% TVAPFParams.M = fs/40;
% TVAPFParams.f_m = 100;
% TVAPFParams.f_b = fs/16;

[yTVAPFMS, yTVAPFWBMat2, TVAPFParams1] = TVAPFSynthesis(fVecCP, env, TVAPFParams, 0, fs);

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
% TVAPFParams.M = 100;
% TVAPFParams.f_m = 250;
% TVAPFParams.f_b = 900;

TVAPFParams.M = 1500;
TVAPFParams.f_m = 100;
TVAPFParams.f_b = 2750;

% loopback FM
BVec = linspace(0.99, 0.989, N);
[yLBFM, ~] = feedbackFMSynthesis(fcVecCP, BVec, env, fs);

% time-varying allpass filter each loopback FM oscillator
yTVAPFLB = zeros(1, N);
envOnes = ones(size(env));
[~, yLBFMMat] = feedbackFMSynthesis(fcVecCP, BVec, envOnes, fs);
for i=1:Nf
    f_pi = fVecCP(i);
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
