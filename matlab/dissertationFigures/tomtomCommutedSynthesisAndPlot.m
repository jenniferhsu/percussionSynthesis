% tomtomCommutedSynthesisAndPlot.m
%
% tomtomSynthesisExample.m should be called to save the tom tom
% resonator files. This file plots and performs commuted synthesis for the 
% tom tom resonator (taiko) created using:
%
%   a. modal synthesis
%   b. the time-varying allpass filtering for sinusoids used in 
%       traditional modal synthesis
%   c. loopback FM modal synthesis
%   d. time-varying allpass filtering on loopback FM oscillators
%
% The excitation is a 0.05 second long noise burst bandpass filtered from
% 120Hz to 4000Hz.  The acoustic resonator IR is a taiko being hit and
% retrieved from freesound.org.

addpath(genpath('../proofOfConcept'));


%% input parameters

fs = 44100;
dur = 2;
N = dur*fs;
savePlots = 1;
writeAudioFiles = 1;

plotSaveDir = 'figures/';
if ~exist(plotSaveDir, 'dir')
    mkdir(plotSaveDir);
end

audioSaveDir = 'audioExamples/tomtom/commutedSynth/';
if ~exist(audioSaveDir, 'dir')
    mkdir(audioSaveDir);
end

yMSWav = 'audioExamples/tomtom/modalSynthesis.wav';
yTVAPFMSWav = 'audioExamples/tomtom/TVAPFModalSynthesis.wav';
yLBWav = 'audioExamples/tomtom/loopbackFM.wav';
yTVAPFLBWav = 'audioExamples/tomtom/TVAPFLoopbackFM.wav';
resIRWav = '../proofOfConcept/resonatorIRs/taiko/taiko2cut.wav';

%% read in audio files

[yMS, ~] = audioread(yMSWav);
[yTVAPFMS, ~] = audioread(yTVAPFMSWav);
[yLBFM, ~] = audioread(yLBWav);
[yTVAPFLB, ~] = audioread(yTVAPFLBWav);
[resIR, ~] = audioread(resIRWav);


%%  generate noise burst excitation

% parameters
durNB = 0.05;
lowFreq = 120;
highFreq = 4000;

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
    
    synthExample = 'Tom tom';

    figure
    subplot(221)
    spectrogram(yMS_NB, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('MS', synthExample));
    ylim([0 5]) 
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(222)
    spectrogram(yLB_NB, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Loopback FM MS', synthExample));
    ylim([0 5]) 
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(223)
    spectrogram(yTVAPFMS_NB, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Time-varying APF MS', synthExample));
    ylim([0 5]) 
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(224)
    spectrogram(yTVAPFLB_NB, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Time-varying APF of Loopback FM MS', synthExample));
    ylim([0 5]) 
    colorbar('off');
    set(gca, 'FontSize', 15);

    sgtitle(sprintf('%s synthesis', synthExample), 'fontsize', 15);

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 6];
    print([plotSaveDir synthExample '_noiseBurst'], '-depsc', '-r0')

end
