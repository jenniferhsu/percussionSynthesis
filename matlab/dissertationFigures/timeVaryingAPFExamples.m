% timeVaryingAPFExamples.m
%
% This creates sound examples to illustrate what the time-varying APF 
% sounds like

addpath(genpath('../proofOfConcept/'));

%% parameters

% input parameters
fs = 44100;
dur = 1;

f0 = 304;
B = 0.9;
g = 0.9999;

saveFigures = 1;
writeAudioFiles = 1;

plotOutDir = 'figures/';
audioOutDir = 'audioExamples/presentations/meeting-presentation/';

envMs = 1; % ms

% derived paramters
N = dur*fs;
T = 1/fs;
nVec = 0:1:N-1;

% parameters
TVAPFParams.M = 5000;
TVAPFParams.f_m = 500;
TVAPFParams.f_b = 1200;

%% set up envelope so that sounds don't click at beginning and end
env = ones(1, N);
envSamps = floor(fs * (envMs/1000));

env(1:envSamps) = linspace(0, 1, envSamps);
env(end-envSamps+1:end) = linspace(1, 0, envSamps);


%% sinusoid

ySine = exp(j*w0*nVec*T);
ySine = ySine .* env;


%% time-varying APF

[yTVAPF, ~, TVAPFParams1] = TVAPFSynthesis(f0, ones(1, N), TVAPFParams, 0, fs);
yTVAPF = yTVAPF .* env;


%% write audio files

if writeAudioFiles
    audiowrite([audioOutDir 'timeVaryingAPFExamplesSinusoid.wav'], scaleForSavingAudio(real(ySine)), fs);
    audiowrite([audioOutDir 'timeVaryingAPFExamplestimeVaryingAPF.wav'], scaleForSavingAudio(real(yTVAPF)), fs);
end


%% plot spectrograms

if saveFigures

    figure
    spectrogram(real(ySine), hann(256), 128, 1024, fs, 'yaxis');
    ylim([0 9])
    colorbar('off')
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([plotOutDir 'timeVaryingAPFExamplesSinusoid'], '-depsc', '-r0')

    figure
    spectrogram(real(yTVAPF), hann(256), 128, 1024, fs, 'yaxis');
    ylim([0 9])
    colorbar('off')
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([plotOutDir 'timeVaryingAPFExamplestimeVaryingAPF'], '-depsc', '-r0')



end