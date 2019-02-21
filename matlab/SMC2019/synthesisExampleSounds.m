% synthesisExampleSounds.m
% this script performs commuted synthesis for the synthesis examples on the
% accompanying website for our SMC2019 paper
% for the audio files for the synthesis examples included in the paper,
% look at synthesisExamplePlots

addpath(genpath('../proofOfConcept'));

%% input parameters
fs = 44100;
N = 2*fs;

%synthExample = 'marimba';
synthExample = 'tom tom';
%synthExample = 'circular plate';
%synthExample = 'bass drum';
%synthExample = 'kettledrum';
%synthExample = 'o-daiko';
%synthExample = 'wood block';

% Marimba example
if strcmp('marimba', synthExample) == 1
    yLBFMWavs = {'../proofOfConcept/audioExamples/marimba/yFBFMMemb4.wav'};
    yMSWav = '../proofOfConcept/audioExamples/marimba/yMS.wav';
    resIRWav = '../proofOfConcept/resonatorIRs/marimbaTube.wav';
end

% Tom tom examples
if strcmp('tom tom', synthExample) == 1
    yLBFMWavs = {'../proofOfConcept/audioExamples/tomtom/yFBFMMemb3.wav', ...
                 '../proofOfConcept/audioExamples/tomtom/ySAPFMemb3.wav'};
                %'../proofOfConcept/audioExamples/tomtom/ySAPFMemb4.wav'
    yMSWav = '../proofOfConcept/audioExamples/tomtom/yMS.wav';
    resIRWav = '../proofOfConcept/resonatorIRs/taiko/taiko2cut.wav';
end

% Circular plate examples
if strcmp('circular plate', synthExample) == 1
    yLBFMWavs = {'../proofOfConcept/audioExamples/circularPlate/yFBFMCP1.wav', ...
                 '../proofOfConcept/audioExamples/circularPlate/yFBFMCP2.wav', ...
                 '../proofOfConcept/audioExamples/circularPlate/yFBFMCP4.wav', ...
                 '../proofOfConcept/audioExamples/circularPlate/ySAPFCP3.wav'};
    resIRWav = 'resonatorIRs/CarpenterCenter.wav';
    yMSWav = '../proofOfConcept/audioExamples/circularPlate/yMS.wav';
    %resIRWav = 'resonatorIRs/3000CStreetGarageStairwell.wav';
end

% Bass drum examples
if strcmp('bass drum', synthExample) == 1
    yLBFMWavs = {'../proofOfConcept/audioExamples/membrane/bass-drum/yFBFMMemb3.wav', ...
                 '../proofOfConcept/audioExamples/membrane/bass-drum/ySAPFMemb3.wav'};
    yMSWav = '../proofOfConcept/audioExamples/membrane/bass-drum/yMS.wav';
    resIRWav = '../proofOfConcept/resonatorIRs/beduk/beduk2cut.wav';
end

% Kettledrum examples
if strcmp('kettledrum', synthExample) == 1
    yLBFMWavs = {'../proofOfConcept/audioExamples/membrane/kettledrum/yFBFMMemb3.wav', ...
                 '../proofOfConcept/audioExamples/membrane/kettledrum/ySAPFMemb3.wav'};
    yMSWav = '../proofOfConcept/audioExamples/membrane/kettledrum/yMS.wav';
    resIRWav = '../proofOfConcept/resonatorIRs/beduk/beduk2cut.wav';
end

% o-daiko examples
if strcmp('o-daiko', synthExample) == 1
    yLBFMWavs = {'../proofOfConcept/audioExamples/membrane/o-daiko/yFBFMMemb3.wav', ...
                 '../proofOfConcept/audioExamples/membrane/o-daiko/ySAPFMemb4.wav'};
    yMSWav = '../proofOfConcept/audioExamples/membrane/o-daiko/yMS.wav';
    resIRWav = '../proofOfConcept/resonatorIRs/beduk/beduk2cut.wav';
end

% Wood block
if strcmp('wood block', synthExample) == 1
    yLBFMWav = '../proofOfConcept/audioExamples/woodBlocks/ySAPFWB2.wav';
    yMSWav = '../proofOfConcept/audioExamples/woodBlocks/yMS.wav';
    resIRWav = '../proofOfConcept/resonatorIRs/113620__vidsyn__miscsoftnaturalgtrloud2-2.wav';
    %resIRWav = '../proofOfConcept/resonatorIRs/taiko/taiko2cut.wav';
end

[yMS, ~] = audioread(yMSWav);

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


%% create excitations and perform commuted synthesis

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

audioOutDir = 'audioExamples/synthesisExampleSounds/';
if ~exist(audioOutDir)
    mkdir(audioOutDir)
end

nFiles = size(yLBFMWavs, 2);
yLBFM_RC = zeros(nFiles,N);
yLBFM_NB = zeros(nFiles,N);

for i=1:nFiles
    
    % raised cosine
    rc = percSynth(dexcitations(:,1), yLBFMWavs{i}, resIRWav);
    yLBFM_RC(i,:) = rc(1:N);

    % noise burst
    nb = percSynth(dexcitations(:,2), yLBFMWavs{i}, resIRWav);
    yLBFM_NB(i,:) = nb(1:N);
    
    % traditional modal synthesis with commuted synthesis
    if i==1
        rc_ms = percSynth(dexcitations(:,1), yMSWav, resIRWav);
        nb_ms = percSynth(dexcitations(:,2), yMSWav, resIRWav);
        audiowrite([audioOutDir synthExample '_MS_RC.wav'], scaleForSavingAudio(rc_ms(1:N)), fs);
        audiowrite([audioOutDir synthExample '_MS_NB.wav'], scaleForSavingAudio(nb_ms(1:N)), fs);
    end
    
    % write audio
    s = strsplit(yLBFMWavs{i}, '/');
    ss = strsplit(s{end}, '.');
    wavname = ss{1};
    
    if strcmp('marimba', synthExample) || ...
            strcmp('circular plate', synthExample) || ... 
            strcmp('bass drum', synthExample) || ... 
            strcmp('kettledrum', synthExample) || ... 
            strcmp('o-daiko', synthExample) || ...
            strcmp('tom tom', synthExample)
        audiowrite([audioOutDir synthExample '_' wavname '_LBFM_RC.wav'], scaleForSavingAudio(yLBFM_RC(i,:)), fs);
    end

    if strcmp('wood block', synthExample) || ...
            strcmp('circular plate', synthExample) || ... 
            strcmp('tom tom', synthExample) || ... 
            strcmp('bass drum', synthExample) || ... 
            strcmp('kettledrum', synthExample) || ... 
            strcmp('o-daiko', synthExample) || ...
            strcmp('tom tom', synthExample) || ...
            strcmp('marimba', synthExample)
        audiowrite([audioOutDir synthExample '_' wavname '_LBFM_NB.wav'], scaleForSavingAudio(yLBFM_NB(i,:)), fs);
    end

end

