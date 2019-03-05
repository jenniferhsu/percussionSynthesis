% modalSynthesisDiagramFigures.m
%
% This script creates the modal synthesis figures (individual sinusoids,
% loopback FM oscillators, envelopes, and outputs)
% that I use for the modal synthesis diagram in my dissertation
% and presentations. Here, we use the modal frequencies from a marimba.

addpath(genpath('../proofOfConcept/'));

%% input parameters

% general
fs = 44100;
dur = 1;

saveFigures = 1;
writeAudioFiles = 1;
plotOutDir = 'figures/'
audioOutDir = 'audioExamples/presentations/meeting-presentation/';

if ~exist(plotOutDir)
    mkdir(plotOutDir)
end

if ~exist(audioOutDir, 'dir')
    mkdir(audioOutDir)
end

% derived parameters

N = fs*dur;
T = 1/fs;

% feedback FM
B = 0.3;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

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

%% modal synthesis with exponential decay using the modes
yMS = zeros(1, N);
yMSMat = zeros(Nf, N);
nVec = 0:T:(dur-T);

for i=1:Nf
    f = fVecBar(i);
    yMSMat(i,:) = exp(1j*2*pi*f*nVec);
    yMS = yMS + (yMSMat(i,:) .* env(i,:));
end


%% loopback FM modal synthesis

% loopback FM
%BVec2 = linspace(0.5, 0.7, N);
BVec2 = linspace(0.5, 0.55, N);
fcVecBar2 = fVecBar./(sqrt(1-BVec2(1)^2));
[yLBFM, yLBFMMat] = feedbackFMSynthesis(fcVecBar2, BVec2, env, fs);

%% time-varying allpass modal synthesis

% parameters
TVAPFParams.M = 100;
TVAPFParams.f_m = 250;
TVAPFParams.f_b = 900;

[yTVAPFMS, yTVAPFMSMat, TVAPFParams1] = TVAPFSynthesis(fVecBar, env, TVAPFParams, 0, fs);


%% time-varying allpass loopback FM

% parameters
TVAPFParams.M = 200;
TVAPFParams.f_m = 300;
TVAPFParams.f_b = 1200;

% time-varying allpass filter each loopback FM oscillator
yTVAPFLB = zeros(1, N);
yTVAPFLBMat = zeros(Nf, N);
envOnes = ones(size(env));
[~, yLBFMMat2] = feedbackFMSynthesis(fcVecBar2, BVec2, envOnes, fs);
for i=1:Nf
    f_pi = fVecBar(i);
    %yTVAPFLB = yTVAPFLB + TVAPF2(yLBFMMat(i,:), f_pi, TVAPFParams1.f_b(i), TVAPFParams1.M(i), TVAPFParams1.f_m(i), fs);
    y1 = TVAPF2(yLBFMMat2(i,:), f_pi, TVAPFParams.f_b, TVAPFParams.M, TVAPFParams.f_m, fs);
    yTVAPFLBMat(i,:) = y1;
    yTVAPFLB = yTVAPFLB + y1 .* env(i,:);
end


%% save audio files

if writeAudioFiles
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% MODAL SYNTHESIS %%% 
    % output
    audiowrite([audioOutDir 'modalSynthesisOutput.wav'], scaleForSavingAudio(real(yMS)), fs);
    
    % individual sinusoids
    for i=1:Nf
        y = real(yMSMat(i,:) .* env(i,:));
        audiowrite([audioOutDir 'modalSynthesisSine' num2str(i) '.wav'], scaleForSavingAudio(y), fs);
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% LOOPBACK FM MS %%% 
    % output
    audiowrite([audioOutDir 'loopbackFMMSOutput.wav'], scaleForSavingAudio(real(yLBFM)), fs);
    
    % loopback FM oscillators
    for i=1:Nf
        y = real(yLBFMMat(i,:) .* env(i,:));
        audiowrite([audioOutDir 'loopbackFMMSOsc' num2str(i) '.wav'], scaleForSavingAudio(y), fs);
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TIME-VARYING APF MS OUTPUT %%%
    audiowrite([audioOutDir 'timeVaryingAPFMSOutput.wav'], scaleForSavingAudio(real(yTVAPFMS)), fs);
    
    % time-varying APF MS oscillators
    for i=1:Nf
        y = real(yTVAPFMSMat(i,:) .* env(i,:));
        audiowrite([audioOutDir 'timeVaryingMSOsc' num2str(i) '.wav'], scaleForSavingAudio(y), fs);
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TIME-VARYING APF LOOPBACK FM OUTPUT %%%
    audiowrite([audioOutDir 'timeVaryingAPFLoopbackFMMSOutput.wav'], scaleForSavingAudio(real(yTVAPFLB)), fs);
    
    % time-varying APF loopback FM MS oscillators
    for i=1:Nf
        y = real(yTVAPFLBMat(i,:) .* env(i,:));
        audiowrite([audioOutDir 'timeVaryingLoopbackFMOsc' num2str(i) '.wav'], scaleForSavingAudio(y), fs);
    end 
    
end


%% plots

if saveFigures
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% MODAL SYNTHESIS %%% 
    
    % sinusoids
    M = 100; % end index
    figure

    for i=1:Nf
        subplot(7, 1, i)
        plot(nVec(1:M), real(yMSMat(i,1:M)), 'linewidth', 2);
        xlim([nVec(1) nVec(M)]);
        set(gca,'linewidth', 3)
        set(gca,'XTick',[], 'YTick', [])
    end

    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    print([plotOutDir 'modalSynthesisDiagramSinusoids'], '-depsc', '-r0')
    
    % MS output
    figure
    plot(nVec, real(yMS), 'linewidth', 2);
    xlim([nVec(1) nVec(end)]);
    ylim([-1 1]);
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([plotOutDir 'modalSynthesisDiagramMSOutput'], '-depsc', '-r0')
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% LOOPBACK FM MS %%% 
    
    % loopback FM oscillators
    M = 100; % end index
    figure

    for i=1:Nf
        subplot(7, 1, i)
        plot(nVec(1:M), real(yLBFMMat(i,1:M)), 'linewidth', 2);
        xlim([nVec(1) nVec(M)]);
        set(gca,'linewidth', 3)
        set(gca,'XTick',[], 'YTick', [])
    end

    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    print([plotOutDir 'modalSynthesisDiagramLoopbackFM'], '-depsc', '-r0')
    
    % Loopback FM output
    figure
    plot(nVec, real(yLBFM), 'linewidth', 2);
    xlim([nVec(1) nVec(end)]);
    ylim([-1 1]);
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([plotOutDir 'modalSynthesisDiagramLoopbackFMOutput'], '-depsc', '-r0')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TIME-VARYING APF MS OUTPUT %%%
    
    % time-varying APF oscillators
    M = 100; % end index
    figure

    for i=1:Nf
        subplot(7, 1, i)
        plot(nVec(1:M), real(yTVAPFMSMat(i,1:M)), 'linewidth', 2);
        xlim([nVec(1) nVec(M)]);
        set(gca,'linewidth', 3)
        set(gca,'XTick',[], 'YTick', [])
    end

    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    print([plotOutDir 'modalSynthesisDiagramTimeVaryingAPF'], '-depsc', '-r0')
    
    % time-varying APF output
    figure
    plot(nVec, real(yTVAPFMS), 'linewidth', 2);
    xlim([nVec(1) nVec(end)]);
    ylim([-1 1]);
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([plotOutDir 'modalSynthesisDiagramTimeVaryingAPFOutput'], '-depsc', '-r0')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TIME-VARYING APF LOOPBACK FM OUTPUT %%%
    
    % time-varying APF loopback FM oscillators
    M = 100; % end index
    figure

    for i=1:Nf
        subplot(7, 1, i)
        plot(nVec(1:M), real(yTVAPFLBMat(i,1:M)), 'linewidth', 2);
        xlim([nVec(1) nVec(M)]);
        set(gca,'linewidth', 3)
        set(gca,'XTick',[], 'YTick', [])
    end

    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    print([plotOutDir 'modalSynthesisDiagramTimeVaryingAPFLoopbackFM'], '-depsc', '-r0')
    
    % time-varying APF loopback FM output
    figure
    plot(nVec, real(yTVAPFLB), 'linewidth', 2);
    xlim([nVec(1) nVec(end)]);
    ylim([-1 1]);
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([plotOutDir 'modalSynthesisDiagramTimeVaryingAPFLoopbackFMOutput'], '-depsc', '-r0')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% EXPONENTIALLY DECAYING AMPLITUDE ENVELOPES %%%
    for i=1:Nf
        subplot(7, 1, i)
        plot(nVec, env(i,:), 'linewidth', 2);
        xlim([nVec(1) nVec(end)]);
        ylim([0 1]);
        set(gca,'linewidth', 3)
        set(gca,'XTick',[], 'YTick', [])
    end
    
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    print([plotOutDir 'modalSynthesisDiagramEnvelopes'], '-depsc', '-r0')
    
    
    
    
end
