% loopbackExamples.m
%
% This creates sound examples to illustrate what loopback FM sounds like

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

w0 = 2*pi*f0;
wc = w0/sqrt(1 - B^2);
BVec = g.^(0:N-1);

%% set up envelope so that sounds don't click at beginning and end
env = ones(1, N);
envSamps = floor(fs * (envMs/1000));

env(1:envSamps) = linspace(0, 1, envSamps);
env(end-envSamps+1:end) = linspace(1, 0, envSamps);


%% sinusoid

ySine = exp(j*w0*nVec*T);
ySine = ySine .* env;


%% static loopback FM

yStaticLB = zeros(1, N);
yStaticLB(1) = 1;

for n=2:N
    yStaticLB(n) = exp(j*wc*T*(1 + B*real(yStaticLB(n-1)))) * yStaticLB(n-1);
end
yStaticLB = yStaticLB .* env;


%% time-varying loopback FM

yTVLB = zeros(1, N);
yTVLB(1) = 1;

for n=2:N
    yTVLB(n) = exp(j*wc*T*(1 + BVec(n)*real(yTVLB(n-1)))) * yTVLB(n-1);
end
yTVLB = yTVLB .* env;

%% write audio files

if writeAudioFiles
    audiowrite([audioOutDir 'loopbackFMExamplesSinusoid.wav'], scaleForSavingAudio(real(ySine)), fs);
    audiowrite([audioOutDir 'loopbackFMExamplesStaticLoopbackFM.wav'], scaleForSavingAudio(yStaticLB), fs);
    audiowrite([audioOutDir 'loopbackFMExamplesTimeVaryingLoopbackFM.wav'], scaleForSavingAudio(yTVLB), fs);
end

%% plot spectrograms

if saveFigures

    figure
    spectrogram(real(ySine), hann(256), 128, 1024, fs, 'yaxis');
    ylim([0 7])
    colorbar('off')
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([plotOutDir 'loopbackFMExamplesSinusoid'], '-depsc', '-r0')

    figure
    spectrogram(real(yStaticLB), hann(256), 128, 1024, fs, 'yaxis');
    ylim([0 7])
    colorbar('off')
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([plotOutDir 'loopbackFMExamplesStaticLoopbackFM'], '-depsc', '-r0')


    figure
    spectrogram(real(yTVLB), hann(256), 128, 1024, fs, 'yaxis');
    ylim([0 7])
    colorbar('off')
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([plotOutDir 'loopbackFMExamplesTimeVaryingLoopbackFM'], '-depsc', '-r0')

end