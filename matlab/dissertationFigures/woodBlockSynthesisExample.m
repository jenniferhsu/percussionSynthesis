% woodBlockSynthesisExample.m
%
% Using the wood block modal frequencies, this script performs:
%   a. modal synthesis
%   b. the time-varying allpass filtering for sinusoids used in 
%       traditional modal synthesis
%   c. loopback FM modal synthesis
%   d. time-varying allpass filtering on loopback FM oscillators
% The wood block modal synthesis signals are saved to
%   timeVaryingAPF/audioExamples/woodBlock directory. The file
%   woodBlockCommutedSynthesisAndPlot.m can be called to perform 
%   commuted synthesis on these saved audio files and to plot the file
%   results.

addpath(genpath('../proofOfConcept'));

%% input parameters

% general
fs = 44100;
dur = 0.5;
plotSpectrograms = 1;
writeAudioFiles = 1;
outDir = 'audioExamples/woodBlock/';

if ~exist(outDir)
    mkdir(outDir)
end

%% derived parameters

% feedback FM
g = 0.9999; % pitch glide coefficient
gWood = 0.9993;

N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient


%% analyze a wood block recording to get the modal frequencies
[wb, fs2] = audioread('../proofOfConcept/audioExamples/woodBlocks/218461__thomasjaunism__wood-block.wav');
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

plot(faxis, WBPosdB);
hold on
for i=1:length(spks)
    plot(faxis(slocs(i)), spks(i), 'r*');
end

fVecWB = faxis(slocs);
fVecWB = sort(fVecWB, 'ascend');
fVecWB = fVecWB(2:end); % drop the lowest one at like 5Hz

Nf = length(fVecWB);


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
eWood = gWood.^(linspace(0, N, N)); % exponential decay for woodblock

envWood = zeros(Nf, N);
for i=1:Nf
    envWood(i,:) = aStart(i) * eWood;
end


%% modal synthesis with exponential decay using the modes
yMS = zeros(1, N);
nVec = 0:T:(dur-T);

for i=1:Nf
    f = fVecWB(i);
    yMS = yMS + (exp(1j*2*pi*f*nVec) .* envWood(i,:));
end


%% time-varying allpass traditional modal synthesis

% parameters
TVAPFParams.M = 1100;
TVAPFParams.f_m = 100;
TVAPFParams.f_b = 2800;

% TVAPFParams.M = fs/40;
% TVAPFParams.f_m = 100;
% TVAPFParams.f_b = fs/16;

[yTVAPFMS, yTVAPFWBMat2, TVAPFParams1] = TVAPFSynthesis(fVecWB, envWood, TVAPFParams, 0, fs);

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
TVAPFParams.M = 1100;
TVAPFParams.f_m = 100;
TVAPFParams.f_b = 2800;

B = 0.9;
b0 = (sqrt(1-B^2) - 1)/B;
fVecWB_w0 = fVecWB*sqrt(1 - B^2);

% stretched APF
[yLBFM, ~] = stretchedAPFSynthesis(fVecWB, b0, envWood, fs, 4*fVecWB, 'linear');
%[yLBFM, yLBFMMat] = stretchedAPFSynthesis(fVecWB, b0, envWood, fs, 2*fVecWB, 'linear');

% time-varying allpass filter each loopback FM oscillator
yTVAPFLB = zeros(1, N);
envOnes = ones(size(envWood));
[~, yLBFMMat] = stretchedAPFSynthesis(fVecWB, b0, envOnes, fs, 4*fVecWB, 'linear');
for i=1:Nf
    f_pi = fVecWB(i);
    %yTVAPFLB = yTVAPFLB + TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams1.f_b(i), TVAPFParams1.M(i), TVAPFParams1.f_m(i), fs);
    y1 = TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams.f_b, TVAPFParams.M, TVAPFParams.f_m, fs);
    yTVAPFLB = yTVAPFLB + y1 .* envWood(i,:);
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
