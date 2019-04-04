% kettledrumSynthesisExamples.m
%
% makes synthesis examples for the o-daiko
addpath(genpath('../proofOfConcept'));

%% input parameters

% general
fs = 44100;
dur = 2;
plotSpectrograms = 1;
writeAudioFiles = 1;
outDir = 'audioExamples/kettledrum/V2/';

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


%% kettledrum modal frequencies
% from Science of Percussion Instruments

% KETTLEDRUM (page 8)
f1 = 150; % fundamental frequency of a kettledrum
fVecMembrane = f1 * [0.63, 1, 1.34, 1.44, 1.66, 1.83, 1.98, 2.20, 2.26, 2.29, 2.55, 2.61, 2.66, 2.89];

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

%% find the pitch glide that I am using in this example
% set the pitch glide with feedback parameter B
b0 = (sqrt(1-B^2) - 1)/B;
fVecMembraneEnd = fVecMembrane - 20*ones(1, Nf);

B1 = -2*b0/(b0^2 + 1);
fc = fVecMembrane(1)/(sqrt(1 - B1^2));
B2 = sqrt(1 - (fVecMembraneEnd(1)/fc)^2);

k = (B2 - B1)/(N - 1);
l = B1 - k;
n = 0:N-1;

BVec = k*n + l;
f0VecMembranePG = fcVecMembrane(:) * sqrt(1 - BVec.^2);

Theta = (fcVecMembrane(:)*T/(2*k)) .* ((k*n + l) .* sqrt(1 - (k*n + l).^2) + asin(k*n + l));


%% loopback FM

[ySAPFMemb3, ~] = stretchedAPFSynthesis(fVecMembrane, b0, env, fs, fVecMembraneEnd, 'linearB');


%% modal synthesis (sinusoidal) with exponential decay using the modes

[yMSPG, ~] = modalSynthSine(Theta, env, fs, 1);


%% modal synthesis (bandpass) with exponential decay using the modes
[yMSBP, ~] = modalSynthBP(f0VecMembranePG, env, fs);

if plotSpectrograms == 1
    figure
    subplot(311)
    spectrogram(real(yMSPG), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, dur, N), f0VecMembranePG/1000, 'r');
    ylim([0 15])
    title('MS spectrogram (sinusoid)')
    subplot(312)
    spectrogram(real(yMSBP), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, dur, N), f0VecMembranePG/1000, 'r');
    ylim([0 15])
    title('MS spectrogram (bandpass)')
    subplot(313)
    spectrogram(real(ySAPFMemb3), hann(256), 128, 1024, fs, 'yaxis');
    hold on
    plot(linspace(0, dur, N), f0VecMembranePG/1000, 'r');
    ylim([0 15])
    title('loopback FM MS spectrogram')
end


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
    audiowrite([outDir 'ySAPFMemb3MSPG.wav'], scaleForSavingAudio(real(yMSPG)), fs);
    %audiowrite([outDir 'TVAPFModalSynthesis.wav'], scaleForSavingAudio(real(yTVAPFMS)), fs);
    
    audiowrite([outDir 'ySAPFMemb3.wav'], scaleForSavingAudio(real(ySAPFMemb3)), fs);
    %audiowrite([outDir 'TVAPFLoopbackFM.wav'], scaleForSavingAudio(real(yTVAPFLB)), fs);
    
    

end
