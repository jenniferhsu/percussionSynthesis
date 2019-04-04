% marimbaSynthesisExampleV2.m
%
% Using the marimba modal frequencies, this script performs:
%   a. modal synthesis as bandpass filters and with pitch glides
%   b. the time-varying allpass filtering for sinusoids used in 
%       traditional modal synthesis
%   c. loopback FM modal synthesis
%   d. time-varying allpass filtering on loopback FM oscillators
% The marimba modal synthesis signals are saved to
%   timeVaryingAPF/audioExamples/marimba directory. The file
%   marimbaCommutedSynthesisAndPlot.m can be called to perform 
%   commuted synthesis on these saved audio files and to plot the file
%   results.

addpath(genpath('../proofOfConcept'));

%% input parameters

% general
fs = 44100;
dur = 2;
plotSpectrograms = 1;
writeAudioFiles = 1;
outDir = 'audioExamples/marimba/V2/';

if ~exist(outDir)
    mkdir(outDir)
end

%% derived parameters

% feedback FM
B = 0.3;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

N = fs*dur;
T = 1/fs;
n = linspace(0, N-1, N);
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient


%% ideal bar modal frequencies
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

e = 0.9999.^(linspace(0, N, N));   % exponential decay

env = zeros(Nf, N);
for i=1:Nf
    env(i,:) = aStart(i) * e;
end

%% modal synthesis (no pitch glide)

[yMS, yMSMat] = modalSynthSine(fVecBar', env, fs, 0);


%% find the pitch glide that I am using in this example
% set the pitch glide with feedback parameter B

f0VecBarPG = fcVecBar(:) .* sqrt(1 - BVec'.^2);


%% loopback FM

% synthesize the loopback FM signal
[yLBFM, ~] = feedbackFMSynthesis(fcVecBar, BVec, env, fs);


%% modal synthesis (sinusoidal) with exponential decay using the modes
u = sqrt(1 - BVec.^2)';
C = 0;                      % constant of integration
Theta = (fcVecBar(:)*T/log(g)) .* (u - atanh(u) + C);
[yMSPG, yMSMat] = modalSynthSine(Theta, env, fs, 1);


%% modal synthesis (bandpass) with exponential decay using the modes
[yMSBP, yMSBPMat] = modalSynthBP(f0VecBarPG, env, fs);


if plotSpectrograms == 1
    figure
    subplot(311)
    spectrogram(real(yMSPG), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, 2, N), f0VecBarPG/1000, 'r');
    ylim([0 3])
    title('MS spectrogram (sinusoid)')
    subplot(312)
    spectrogram(real(yMSBP), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, 2, N), f0VecBarPG/1000, 'r');
    ylim([0 3])
    title('MS spectrogram (bandpass)')
    subplot(313)
    spectrogram(real(yLBFM), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, 2, N), f0VecBarPG/1000, 'r');
    ylim([0 3])
    title('loopback FM MS spectrogram')
end


%% Other loopback FM examples for the webpage
B1 = 0.7;
B2 = 0.25;
k = (B2 - B1)/(N - 1);
l = B1 - k;
%BVec3 = linspace(0.7, 0.25, N);
BVec3 = k*n + l;

fcVecBar3 = fVecBar./(sqrt(1-BVec3(1)^2));
[yLBFMMemb4, yFBFMMembMat4] = feedbackFMSynthesis(fcVecBar3, BVec3, env, fs);


%% find the pitch glide that I am using in this example
f0VecBar3 = fcVecBar3(:) .* sqrt(1 - BVec3.^2);


%% modal synthesis (sinusoidal) with exponential decay using the modes
f0VecBar3PG = (fcVecBar3(:)*T/(2*k)) .* ((k*n + l) .* sqrt(1 - (k*n + l).^2) + asin(k*n + l));
[yLBFMMemb4MSPG, ~] = modalSynthSine(f0VecBar3PG, env, fs, 1);

%% time-varying allpass traditional modal synthesis
% [for dissertation, not SMC]

% % parameters
% TVAPFParams.M = 100;
% TVAPFParams.f_m = 250;
% TVAPFParams.f_b = 900;
% 
% [yTVAPFMS, yTVAPFWBMat2, TVAPFParams1] = TVAPFSynthesis(fVecBar, env, TVAPFParams, 0, fs);
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



%% time-varying allpass loopback FM
% [for dissertation, not SMC]

% % parameters
% TVAPFParams.M = 200;
% TVAPFParams.f_m = 300;
% TVAPFParams.f_b = 1200;
% 
% % time-varying allpass filter each loopback FM oscillator
% yTVAPFLB = zeros(1, N);
% envOnes = ones(size(env));
% [~, yLBFMMat] = feedbackFMSynthesis(fcVecBar2, BVec2, envOnes, fs);
% for i=1:Nf
%     f_pi = fVecBar(i);
%     %yTVAPFLB = yTVAPFLB + TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams1.f_b(i), TVAPFParams1.M(i), TVAPFParams1.f_m(i), fs);
%     y1 = TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams.f_b, TVAPFParams.M, TVAPFParams.f_m, fs);
%     yTVAPFLB = yTVAPFLB + y1 .* env(i,:);
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

%% Linear pitch glide for stretched allpass 
%[for reference]
% BLin = k*n + l;
% ThetaH = (wc*T/(2*k)) * (BLin .* sqrt(1 - BLin.^2) + asin(BLin));
% b0 = (sqrt(1-B.^2) - 1)./B;
% H = (b0 + exp(1j*ThetaH)) ./ (1 + b0.*exp(1j*ThetaH));
% 
% %H = (b0 + exp(1j*cumsum(w0)*T)) ./ (1 + b0.*exp(1j*cumsum(w0)*T));
% % ^^ that's how you do it experimentally
% spectrogram(real(H), hann(256), 128, 1024, fs, 'yaxis');
% hold on
% plot(linspace(0, 2, N), w0/(2*pi*1000), 'r');



%% Write audio files

if writeAudioFiles == 1
    audiowrite([outDir 'modalSynthesis.wav'], scaleForSavingAudio(real(yMS)), fs);
    audiowrite([outDir 'modalSynthesisPitchGlide.wav'], scaleForSavingAudio(real(yMSPG)), fs);
    audiowrite([outDir 'yLBFMMemb4ModalSynthesisPitchGlide.wav'], scaleForSavingAudio(real(yLBFMMemb4MSPG)), fs);
    %audiowrite([outDir 'TVAPFModalSynthesis.wav'], scaleForSavingAudio(real(yTVAPFMS)), fs);
    
    audiowrite([outDir 'loopbackFM.wav'], scaleForSavingAudio(real(yLBFM)), fs);
    audiowrite([outDir 'yLBFMMemb4.wav'], scaleForSavingAudio(real(yLBFMMemb4)), fs);
    %audiowrite([outDir 'TVAPFLoopbackFM.wav'], scaleForSavingAudio(real(yTVAPFLB)), fs);
end
