% circularPlateSynthesisExampleV2.m
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
dur = 2;
plotSpectrograms = 1;
writeAudioFiles = 1;
outDir = 'audioExamples/circularPlate/V2/';

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


%% modal synthesis (no pitch glide)
[yMS, yMSMat] = modalSynthSine(fVecCP', env, fs, 0);



%% loopback FM
%BVec = linspace(0.99, 0.989, N);
n = 0:N-1;
%B1 = 0.99;
%B2 = 0.989;
B1 = 0.91;
B2 = 0.90;
k = (B2 - B1)/(N - 1);
l = B1 - k;
BVec = k*n + l;
[yLBFM, ~] = feedbackFMSynthesis(fcVecCP, BVec, env, fs);


%% find the pitch glide that I am using in my examples
f0VecCP = fcVecCP(:) .* sqrt(1 - BVec.^2);


%% modal synthesis (sinusoidal) with exponential decay using the modes
Thetaf0 = (fcVecCP(:)*T/(2*k)) .* ((k*n + l) .* sqrt(1 - (k*n + l).^2) + asin(k*n + l));
[yMSPG, yMSMat] = modalSynthSine(Thetaf0, env, fs, 1);


%% modal synthesis (bandpass) with exponential decay using the modes
[yMSBP, yMSBPMat] = modalSynthBP(f0VecCP, env, fs);


if plotSpectrograms == 1
    figure
    subplot(311)
    spectrogram(real(yMSPG), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, dur, N), f0VecCP/1000, 'r');
    ylim([0 15])
    title('MS spectrogram (sinusoid)')
    subplot(312)
    spectrogram(real(yMSBP), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, dur, N), f0VecCP/1000, 'r');
    ylim([0 15])
    title('MS spectrogram (bandpass)')
    subplot(313)
    spectrogram(real(yLBFM), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, dur, N), f0VecCP/1000, 'r');
    ylim([0 15])
    title('loopback FM MS spectrogram')
end

%% Other Loopback FM ones for sound examples

% yFBFMCP1
[yFBFMCP1, yFBFMCPMat1] = feedbackFMSynthesis(fVecCP, B, env, fs);

% MS pitch glide (this is a dummy one because it's the same!) - just
% include the one without the pitch glide
f00VecCP = fVecCP * sqrt(1 - B^2);
[yFBFMCP1MSPG, yMSMat] = modalSynthSine(f00VecCP', env, fs, 0);

% yFBFMCP2 
[yFBFMCP2, yFBFMCPMat2] = feedbackFMSynthesis(fcVecCP, B, env, fs);

% MS pitch glide (this is a dummy one because there isn't a pitch glide!) -
% just include the one without the pitch glide, leave the pitch glide MS
% section blank
[yFBFMCP2MSPG, yMSMat] = modalSynthSine(fVecCP', env, fs, 0);

% ySAPFCP3
b0 = (sqrt(1-B^2) - 1)/B;
BVec = linspace(0.99, 0.989, N);
fcVecCPStart = fcVecCP .* sqrt(1 - BVec(1)^2);
fcVecCPEnd = fcVecCP .* sqrt(1 - BVec(end)^2);
[ySAPFCP3, ySAPFCPMat3] = stretchedAPFSynthesis(fcVecCPStart, b0, env, fs, fcVecCPEnd, 'linearB');

% MS pitch glide
B1 = -2*b0/(b0^2 + 1);

fc = fcVecCPStart(1)/(sqrt(1 - B1^2));
B2 = sqrt(1 - (fcVecCPEnd(1)/fc)^2);

k = (B2 - B1)/(N - 1);
l = B1 - k;
n = 0:N-1;

Theta = (2*pi) * (fcVecCPStart(:)*T/(2*k)) .* ((k*n + l) .* sqrt(1 - (k*n + l).^2) + asin(k*n + l));
[ySAPFCP3MSPG, ~] = modalSynthSine(Theta, env, fs, 1);


%% time-varying allpass traditional modal synthesis
% [for dissertation, not SMC]

% parameters
% TVAPFParams.M = 11500;
% TVAPFParams.f_m = 1000;
% TVAPFParams.f_b = 2750;
% 
% [yTVAPFMS, yTVAPFWBMat2, TVAPFParams1] = TVAPFSynthesis(fVecCP, env, TVAPFParams, 0, fs);
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
% 
% % time-varying allpass filter each loopback FM oscillator
% yTVAPFLB = zeros(1, N);
% envOnes = ones(size(env));
% [~, yLBFMMat] = feedbackFMSynthesis(fcVecCP, BVec, envOnes, fs);
% for i=1:Nf
%     f_pi = fVecCP(i);
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
    
    audiowrite([outDir 'loopbackFM.wav'], scaleForSavingAudio(real(yLBFM)), fs);
    %audiowrite([outDir 'TVAPFLoopbackFM.wav'], scaleForSavingAudio(real(yTVAPFLB)), fs);
    
    audiowrite([outDir 'yFBFMCP1.wav'], scaleForSavingAudio(real(yFBFMCP1)), fs);
    audiowrite([outDir 'yFBFMCP1MSPG.wav'], scaleForSavingAudio(real(yFBFMCP1MSPG)), fs);
    
    audiowrite([outDir 'yFBFMCP2.wav'], scaleForSavingAudio(real(yFBFMCP2)), fs);
    audiowrite([outDir 'yFBFMCP2MSPG.wav'], scaleForSavingAudio(real(yFBFMCP2MSPG)), fs);
    
    audiowrite([outDir 'ySAPFCP3.wav'], scaleForSavingAudio(real(ySAPFCP3)), fs);
    audiowrite([outDir 'ySAPFCP3MSPG.wav'], scaleForSavingAudio(real(ySAPFCP3MSPG)), fs);
    

end
