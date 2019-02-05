% synthesisExamplePlots.m
% this script performs commuted synthesis for the synthesis examples in the
% 2019 SMC paper

addpath(genpath('../proofOfConcept'));

%% input parameters
fs = 44100;
N = 2*fs;

%synthExample = 'marimba';
%synthExample = 'tom tom';
%synthExample = 'circular plate';
synthExample = 'wood block';

% Marimba example
if strcmp('marimba', synthExample) == 1
    yLBFMWav = '../proofOfConcept/audioExamples/marimba/yFBFMMemb3.wav';
    yMSWav = '../proofOfConcept/audioExamples/marimba/yMS.wav';
    resIRWav = '../proofOfConcept/resonatorIRs/marimbaTube.wav';
end

% Tom tom examples
if strcmp('tom tom', synthExample) == 1
    yLBFMWav = '../proofOfConcept/audioExamples/tomtom/ySAPFMemb4.wav';
    yMSWav = '../proofOfConcept/audioExamples/tomtom/yMS.wav';
    resIRWav = '../proofOfConcept/resonatorIRs/taiko/taiko2cut.wav';
end

% Circular plate example
if strcmp('circular plate', synthExample) == 1
    yLBFMWav = '../proofOfConcept/audioExamples/circularPlate/yFBFMCP3.wav';
    yMSWav = '../proofOfConcept/audioExamples/circularPlate/yMS.wav';
    resIRWav = 'resonatorIRs/CarpenterCenter.wav';
end

% Wood block
if strcmp('wood block', synthExample) == 1
    yLBFMWav = '../proofOfConcept/audioExamples/woodBlocks/ySAPFWB2.wav';
    yMSWav = '../proofOfConcept/audioExamples/woodBlocks/yMS.wav';
    resIRWav = '../proofOfConcept/resonatorIRs/113620__vidsyn__miscsoftnaturalgtrloud2-2.wav';
    %resIRWav = '../proofOfConcept/resonatorIRs/taiko/taiko2cut.wav';
end

[yLBFM, ~] = audioread(yLBFMWav);
[yMS, ~] = audioread(yMSWav);
[resIR, ~] = audioread(resIRWav);

% modal synthesis signal location 
%yMSWav = '../proofOfConcept/audioExamples/membrane/kettledrum/ySAPFMemb2.wav';
%yMSWav = '../proofOfConcept/audioExamples/membrane/bass-drum/ySAPFMemb2.wav';
%yMSWav = '../proofOfConcept/audioExamples/membrane/o-daiko/ySAPFMemb3.wav';

% resonating body impulse response location
%resIRWav = '../proofOfConcept/resonatorIRs/195790__klankbeeld__cinematic-boom-130730-06.wav';
%resIRWav = '../proofOfConcept/resonatorIRs/10543__batchku__hit-low-14-002.aiff';
%resIRWav = '../proofOfConcept/resonatorIRs/39303__the-semen-incident__pvc-tub-bass.wav';
%resIRWav = '../proofOfConcept/resonatorIRs/45512__tim-kahn__whipper-1.flac';
%resIRWav = '../proofOfConcept/resonatorIRs/45595__tim-kahn__futon-bd-1.aiff';
%resIRWav = '../proofOfConcept/resonatorIRs/45596__tim-kahn__futon-bd-2.aiff';
%resIRWav = '../proofOfConcept/resonatorIRs/snareShell.wav';
%resIRWav = '../proofOfConcept/resonatorIRs/taiko/taiko2.wav';

% Room IRs
%resIRWav = 'resonatorIRs/3000CStreetGarageStairwell.wav';


%% create excitations and perform convolutional synthesis

excitations = zeros(N,2);
dexcitations = zeros(N,2);  % first derivative of excitations

%% generate raised cosine excitation

winLength = 8;

% raised cosine/Hann window
n = winLength/2:winLength-1;
w = 0.5 * (1 - cos((2*pi*n)/(winLength-1)));

excitations(1:winLength/2, 1) = ones(winLength/2, 1);
excitations(winLength/2+1:winLength, 1) = w;

% take derivative for velocity
dexcitations(:,1) = [diff(excitations(:,1)); 0]; 


%% generate filtered noise burst excitation

durNB = 0.05;
lowFreq = 120;
highFreq = 4000;
if strcmp('wood block', synthExample)
    durNB = .001;
    highFreq = 12000;
end

% noise burst
sampNB = ceil(durNB*fs);
excitations(:,2) = [(2*rand(1, sampNB))-1 zeros(1, N-sampNB)]';

[B, A] = butter(5, [lowFreq/(fs/2) highFreq/(fs/2)], 'bandpass');
%freqz(B,A)
excitations(:,2) = filter(B, A, excitations(:,2));

% take derivative for velocity
dexcitations(:,2) = [diff(excitations(:,2)); 0]; 


%% synthesis

% basic modal synthesis
yMS_RC = percSynth(dexcitations(:,1), yMSWav, resIRWav);
yMS_RC = yMS_RC(1:N);

yMS_NB = percSynth(dexcitations(:,2), yMSWav, resIRWav);
yMS_NB = yMS_NB(1:N);

% LBFM synthesis
yLBFM_RC = percSynth(dexcitations(:,1), yLBFMWav, resIRWav);
yLBFM_RC = yLBFM_RC(1:N);

yLBFM_NB = percSynth(dexcitations(:,2), yLBFMWav, resIRWav);
yLBFM_NB = yLBFM_NB(1:N);

%% plots

saveDir = 'figures/';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% raised cosine
figure
subplot(211)
spectrogram(yMS_RC, hann(256), 128, 1024, fs, 'yaxis');
title(sprintf('modal synthesis of %s with raised cosine excitation', synthExample));
set(gca, 'FontSize', 15);
if strcmp('tom tom', synthExample)
    ylim([0 5]) 
else
    ylim([0 18])
end
if strcmp('wood block', synthExample)
    xlim([0 0.5]);
end
subplot(212)
spectrogram(yLBFM_RC, hann(256), 128, 1024, fs, 'yaxis');
title(sprintf('LBFM modal synthesis of %s with raised cosine excitation', synthExample));
set(gca, 'FontSize', 15);
if strcmp('tom tom', synthExample)
    ylim([0 5]) 
else
    ylim([0 18])
end
if strcmp('wood block', synthExample)
    xlim([0 0.5]);
end
saveas(gcf, [saveDir synthExample '_raisedCosine'], 'epsc')


% noise burst
figure
subplot(211)
spectrogram(yMS_NB, hann(256), 128, 1024, fs, 'yaxis');
title(sprintf('modal synthesis of %s with filtered noise burst', synthExample));
set(gca, 'FontSize', 15);
if strcmp('tom tom', synthExample)
    ylim([0 5]) 
else
    ylim([0 18])
end
if strcmp('wood block', synthExample)
    xlim([0 0.5]);
end
subplot(212)
spectrogram(yLBFM_NB, hann(256), 128, 1024, fs, 'yaxis');
title(sprintf('LBFM modal synthesis of %s with filtered noise burst', synthExample));
set(gca, 'FontSize', 15);
if strcmp('tom tom', synthExample)
    ylim([0 5])
else
    ylim([0 18])
end
if strcmp('wood block', synthExample)
    xlim([0 0.5]);
end
saveas(gcf, [saveDir synthExample '_noiseBurst'], 'epsc')

%% write audio
audioOutDir = 'audioExamples/synthesisExamples/';
if ~exist(audioOutDir)
    mkdir(audioOutDir)
end

if strcmp('marimba', synthExample) || ...
        strcmp('circular plate', synthExample)
    audiowrite([audioOutDir synthExample '_MS_RC.wav'], scaleForSavingAudio(yMS_RC), fs);
    audiowrite([audioOutDir synthExample '_LBFM_RC.wav'], scaleForSavingAudio(yLBFM_RC), fs);
end

if strcmp('wood block', synthExample) || ...
        strcmp('tom tom', synthExample)
    audiowrite([audioOutDir synthExample '_MS_NB.wav'], scaleForSavingAudio(yMS_NB), fs);
    audiowrite([audioOutDir synthExample '_LBFM_NB.wav'], scaleForSavingAudio(yLBFM_NB), fs);
end