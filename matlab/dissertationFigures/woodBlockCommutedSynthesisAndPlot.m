% woodBlockCommutedSynthesisAndPlot.m
%
% woodBlockSynthesisExample.m should be called to save the wood block
% resonator files. This file plots and performs commuted synthesis for the 
% wood block resonator created using:
%
%   a. modal synthesis
%   b. the time-varying allpass filtering for sinusoids used in 
%       traditional modal synthesis
%   c. loopback FM modal synthesis
%   d. time-varying allpass filtering on loopback FM oscillators
%
% The excitation is a 0.001 second long noise burst bandpass filtered from
% 120Hz to 12000Hz.  The acoustic resonator IR is an acoustic tube created using traditional modal
% synthesis.

addpath(genpath('../proofOfConcept'));


%% input parameters

fs = 44100;
dur = 0.5;
N = dur*fs;
savePlots = 1;
writeAudioFiles = 1;

plotSaveDir = 'figures/';
if ~exist(plotSaveDir, 'dir')
    mkdir(plotSaveDir);
end

audioSaveDir = 'audioExamples/woodBlock/commutedSynth/';
if ~exist(audioSaveDir, 'dir')
    mkdir(audioSaveDir);
end

yMSWav = 'audioExamples/woodBlock/modalSynthesis.wav';
yTVAPFMSWav = 'audioExamples/woodBlock/TVAPFModalSynthesis.wav';
yLBWav = 'audioExamples/woodBlock/loopbackFM.wav';
yTVAPFLBWav = 'audioExamples/woodBlock/TVAPFLoopbackFM.wav';
resIRWav = '../proofOfConcept/resonatorIRs/113620__vidsyn__miscsoftnaturalgtrloud2-2.wav';

%% read in audio files

[yMS, ~] = audioread(yMSWav);
[yTVAPFMS, ~] = audioread(yTVAPFMSWav);
[yLBFM, ~] = audioread(yLBWav);
[yTVAPFLB, ~] = audioread(yTVAPFLBWav);
[resIR, ~] = audioread(resIRWav);


%%  generate raised cosine excitation

% parameters
durNB = .001;
lowFreq = 120;
highFreq = 12000;

% noise burst
sampNB = ceil(durNB*fs);
excitation = [(2*rand(1, sampNB))-1 zeros(1, N-sampNB)]';

[B, A] = butter(5, [lowFreq/(fs/2) highFreq/(fs/2)], 'bandpass');
%freqz(B,A)
excitation = filter(B, A, excitation);

% take derivative for velocity
dexcitation = [diff(excitation); 0]; 


%% commuted synthesis

% traditional modal synthesis
yMS_NB = percSynth(dexcitation, yMSWav, resIRWav);
yMS_NB = yMS_NB(1:N);

% time-varying APF with traditional modal synthesis
yTVAPFMS_NB = percSynth(dexcitation, yTVAPFMSWav, resIRWav);
yTVAPFMS_NB = yTVAPFMS_NB(1:N);

% LBFM synthesis
yLB_NB = percSynth(dexcitation, yLBWav, resIRWav);
yLB_NB = yLB_NB(1:N);

% LBFM synthesis
yTVAPFLB_NB = percSynth(dexcitation, yTVAPFLBWav, resIRWav);
yTVAPFLB_NB = yTVAPFLB_NB(1:N);


%% save audio examples

if writeAudioFiles == 1
    audiowrite([audioSaveDir 'modalSynthesis.wav'], scaleForSavingAudio(real(yMS_NB)), fs);
    audiowrite([audioSaveDir 'TVAPFModalSynthesis.wav'], scaleForSavingAudio(real(yTVAPFMS_NB)), fs);
    
    audiowrite([audioSaveDir 'loopbackFM.wav'], scaleForSavingAudio(real(yLB_NB)), fs);
    audiowrite([audioSaveDir 'TVAPFLoopbackFM.wav'], scaleForSavingAudio(real(yTVAPFLB_NB)), fs);
end

%% plot

if savePlots == 1
    
    synthExample = 'Wood block';

    figure
    subplot(221)
    spectrogram(yMS_NB, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('MS', synthExample));
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(222)
    spectrogram(yLB_NB, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Loopback FM MS', synthExample));
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(223)
    spectrogram(yTVAPFMS_NB, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Time-varying APF MS', synthExample));
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(224)
    spectrogram(yTVAPFLB_NB, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Time-varying APF of Loopback FM MS', synthExample));
    colorbar('off');
    set(gca, 'FontSize', 15);

    sgtitle(sprintf('%s synthesis', synthExample), 'fontsize', 15);

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 6];
    print([plotSaveDir synthExample '_noiseBurst'], '-depsc', '-r0')

end
