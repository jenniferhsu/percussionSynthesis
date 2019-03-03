% circularPlateCommutedSynthesisAndPlot.m
%
% circularPlateSynthesisExample.m should be called to save the circular 
% plate resonator files. This file plots and performs commuted synthesis 
% for the circular plate resonator (room IR) created using:
%
%   a. modal synthesis
%   b. the time-varying allpass filtering for sinusoids used in 
%       traditional modal synthesis
%   c. loopback FM modal synthesis
%   d. time-varying allpass filtering on loopback FM oscillators
%
% The excitation is an 8-sample long raised cosine.  The acoustic 
% resonator IR is a room impulse response retrieved from EchoThief.com.

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

audioSaveDir = 'audioExamples/circularPlate/commutedSynth/';
if ~exist(audioSaveDir, 'dir')
    mkdir(audioSaveDir);
end

yMSWav = 'audioExamples/circularPlate/modalSynthesis.wav';
yTVAPFMSWav = 'audioExamples/circularPlate/TVAPFModalSynthesis.wav';
yLBWav = 'audioExamples/circularPlate/loopbackFM.wav';
yTVAPFLBWav = 'audioExamples/circularPlate/TVAPFLoopbackFM.wav';
resIRWav = 'resonatorIRs/CarpenterCenter.wav';

%% read in audio files

[yMS, ~] = audioread(yMSWav);
[yTVAPFMS, ~] = audioread(yTVAPFMSWav);
[yLBFM, ~] = audioread(yLBWav);
[yTVAPFLB, ~] = audioread(yTVAPFLBWav);
[resIR, ~] = audioread(resIRWav);


%%  generate raised cosine excitation

winLength = 8;
excitation = zeros(N, 1);

% raised cosine/Hann window
n = winLength/2:winLength-1;
w = 0.5 * (1 - cos((2*pi*n)/(winLength-1)));

excitation(1:winLength/2) = ones(winLength/2, 1);
excitation(winLength/2+1:winLength) = w;

% take derivative for velocity
dexcitation = [diff(excitation); 0]; 


%% commuted synthesis

% traditional modal synthesis
yMS_RC = percSynth(dexcitation, yMSWav, resIRWav);
yMS_RC = yMS_RC(1:N);

% time-varying APF with traditional modal synthesis
yTVAPFMS_RC = percSynth(dexcitation, yTVAPFMSWav, resIRWav);
yTVAPFMS_RC = yTVAPFMS_RC(1:N);

% LBFM synthesis
yLB_RC = percSynth(dexcitation, yLBWav, resIRWav);
yLB_RC = yLB_RC(1:N);

% LBFM synthesis
yTVAPFLB_RC = percSynth(dexcitation, yTVAPFLBWav, resIRWav);
yTVAPFLB_RC = yTVAPFLB_RC(1:N);


%% save audio examples

if writeAudioFiles == 1
    audiowrite([audioSaveDir 'modalSynthesis.wav'], scaleForSavingAudio(real(yMS_RC)), fs);
    audiowrite([audioSaveDir 'TVAPFModalSynthesis.wav'], scaleForSavingAudio(real(yTVAPFMS_RC)), fs);
    
    audiowrite([audioSaveDir 'loopbackFM.wav'], scaleForSavingAudio(real(yLB_RC)), fs);
    audiowrite([audioSaveDir 'TVAPFLoopbackFM.wav'], scaleForSavingAudio(real(yTVAPFLB_RC)), fs);
end

%% plot

if savePlots == 1
    
    synthExample = 'Circular plate';

    figure
    subplot(221)
    spectrogram(yMS_RC, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('MS', synthExample));
    ylim([0 18]) 
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(222)
    spectrogram(yLB_RC, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Loopback FM MS', synthExample));
    ylim([0 18]) 
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(223)
    spectrogram(yTVAPFMS_RC, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Time-varying APF MS', synthExample));
    ylim([0 18]) 
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(224)
    spectrogram(yTVAPFLB_RC, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Time-varying APF of Loopback FM MS', synthExample));
    ylim([0 18]) 
    colorbar('off');
    set(gca, 'FontSize', 15);

    sgtitle(sprintf('%s synthesis', synthExample), 'fontsize', 15);

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 6];
    print([plotSaveDir synthExample '_raisedCosine'], '-depsc', '-r0')

end
