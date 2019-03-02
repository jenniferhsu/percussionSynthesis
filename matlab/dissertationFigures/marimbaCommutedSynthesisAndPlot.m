% marimbaCommutedSynthesisAndPlot.m
%
% marimbaSynthesisExample.m should be called to save the marimba resonator
% files. This file plots and performs commuted synthesis for the marimba 
% resonator created using:
%
%   a. modal synthesis
%   b. the time-varying allpass filtering for sinusoids used in 
%       traditional modal synthesis
%   c. loopback FM modal synthesis
%   d. time-varying allpass filtering on loopback FM oscillators
%
% The excitation is an 8-sample long raised cosine excitation and the
% acoustic resonator IR is an acoustic tube created using traditional modal
% synthesis.

addpath(genpath('../proofOfConcept'));


%% input parameters

fs = 44100;
N = 2*fs;
savePlots = 1;
writeAudioFiles = 1;

plotSaveDir = 'figures/';
if ~exist(plotSaveDir, 'dir')
    mkdir(plotSaveDir);
end

audioSaveDir = 'audioExamples/marimba/commutedSynth/';
if ~exist(audioSaveDir, 'dir')
    mkdir(audioSaveDir);
end

yMSWav = 'audioExamples/marimba/modalSynthesis.wav';
yTVAPFMSWav = 'audioExamples/marimba/TVAPFModalSynthesis.wav';
yLBWav = 'audioExamples/marimba/loopbackFM.wav';
yTVAPFLBWav = 'audioExamples/marimba/TVAPFLoopbackFM.wav';
resIRWav = '../proofOfConcept/resonatorIRs/marimbaTube.wav';

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
    
    synthExample = 'Marimba';

    figure
    subplot(221)
    spectrogram(yMS_RC, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('MS', synthExample));
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(222)
    spectrogram(yLB_RC, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Loopback FM MS', synthExample));
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(223)
    spectrogram(yTVAPFMS_RC, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Time-varying APF MS', synthExample));
    colorbar('off');
    set(gca, 'FontSize', 15);

    subplot(224)
    spectrogram(yTVAPFLB_RC, hann(256), 128, 1024, fs, 'yaxis');
    title(sprintf('Time-varying APF of Loopback FM MS', synthExample));
    colorbar('off');
    set(gca, 'FontSize', 15);

    sgtitle(sprintf('%s synthesis', synthExample), 'fontsize', 15);

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 6];
    print([plotSaveDir synthExample '_raisedCosine'], '-depsc', '-r0')

end
