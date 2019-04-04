% woodBlockSynthesisExampleV2.m
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
outDir = 'audioExamples/woodBlock/V2/';

if ~exist(outDir)
    mkdir(outDir)
end

%% derived parameters

% feedback FM
g = 0.9999; % pitch glide coefficient
gWood = 0.999;

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


%% modal synthesis (no pitch glide)
[yMS, yMSMat] = modalSynthSine(fVecWB', envWood, fs, 0);


%% loopback FM
B = 0.9;
b0 = (sqrt(1-B^2) - 1)/B;
fVecWBEnd = 2*fVecWB;

% stretched APF
[yLBFM, ~] = stretchedAPFSynthesis(fVecWB, b0, envWood, fs, fVecWBEnd, 'linearB');
%[yLBFM, yLBFMMat] = stretchedAPFSynthesis(fVecWB, b0, envWood, fs, 2*fVecWB, 'linear');


%% find the pitch glide that I am using in my examples
B1 = -2*b0/(b0^2 + 1);
fc = fVecWB(1)/(sqrt(1 - B1^2));
B2 = sqrt(1 - (fVecWBEnd(1)/fc)^2);

k = (B2 - B1)/(N - 1);
l = B1 - k;
n = 0:N-1;
BVec = k*n + l;

fcVecWB = fVecWB/(sqrt(1 - B1^2));
f0VecWB = fcVecWB(:) .* sqrt(1 - BVec.^2);  % pitch glide to compare to

%% modal synthesis (sinusoidal) with exponential decay using the modes
f0VecWB2 = (fcVecWB(:)*T/(2*k)) .* ((k*n + l) .* sqrt(1 - (k*n + l).^2) + asin(k*n + l));
[yMSPG, yMSMat] = modalSynthSine(f0VecWB2, envWood, fs, 1);


%% modal synthesis (bandpass) with exponential decay using the modes
[yMSBP, yMSBPMat] = modalSynthBP(f0VecWB, envWood, fs);

if plotSpectrograms == 1
    figure
    subplot(311)
    spectrogram(real(yMSPG), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, 2, N), f0VecWB/1000, 'r');
    ylim([0 3])
    title('MS spectrogram (sinusoid)')
    subplot(312)
    spectrogram(real(yMSBP), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, 2, N), f0VecWB/1000, 'r');
    ylim([0 3])
    title('MS spectrogram (bandpass)')
    subplot(313)
    spectrogram(real(yLBFM), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, 2, N), f0VecWB/1000, 'r');
    ylim([0 3])
    title('loopback FM MS spectrogram')
end



%% time-varying allpass traditional modal synthesis
% [for dissertation, not SMC]

% parameters
% TVAPFParams.M = 1100;
% TVAPFParams.f_m = 100;
% TVAPFParams.f_b = 2800;
% 
% [yTVAPFMS, yTVAPFWBMat2, TVAPFParams1] = TVAPFSynthesis(fVecWB, envWood, TVAPFParams, 0, fs);
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
% TVAPFParams.M = 1100;
% TVAPFParams.f_m = 100;
% TVAPFParams.f_b = 2800;
% 
% B = 0.9;
% b0 = (sqrt(1-B^2) - 1)/B;
% fVecWB_w0 = fVecWB*sqrt(1 - B^2);
% 
% % stretched APF
% [yLBFM, ~] = stretchedAPFSynthesis(fVecWB, b0, envWood, fs, 4*fVecWB, 'linear');
% %[yLBFM, yLBFMMat] = stretchedAPFSynthesis(fVecWB, b0, envWood, fs, 2*fVecWB, 'linear');
% 
% % time-varying allpass filter each loopback FM oscillator
% yTVAPFLB = zeros(1, N);
% envOnes = ones(size(envWood));
% [~, yLBFMMat] = stretchedAPFSynthesis(fVecWB, b0, envOnes, fs, 4*fVecWB, 'linear');
% for i=1:Nf
%     f_pi = fVecWB(i);
%     %yTVAPFLB = yTVAPFLB + TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams1.f_b(i), TVAPFParams1.M(i), TVAPFParams1.f_m(i), fs);
%     y1 = TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams.f_b, TVAPFParams.M, TVAPFParams.f_m, fs);
%     yTVAPFLB = yTVAPFLB + y1 .* envWood(i,:);
% end
% 
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
    
    audiowrite([outDir 'loopbackFM.wav'], scaleForSavingAudio(real(yLBFM)), fs);
    %audiowrite([outDir 'TVAPFLoopbackFM.wav'], scaleForSavingAudio(real(yTVAPFLB)), fs);
end
