% synthesisExampleSoundsV2.m
% this script performs commuted synthesis for the synthesis examples on the
% accompanying website for our SMC2019 paper
% for the audio files for the synthesis examples included in the paper,
% look at synthesisExamplePlots



%% input parameters
fs = 44100;
N = 2*fs;

baseDir = '~/Documents/ucsd/dissertation/percussionSynthesis/matlab/SMC2019/';
%synthExample = 'marimba';
synthExample = 'tom tom';
%synthExample = 'circular plate';
%synthExample = 'bass drum';
%synthExample = 'kettledrum';
%synthExample = 'o-daiko';
%synthExample = 'wood block';

% Marimba example
if strcmp('marimba', synthExample) == 1
    yLBFMWavs = {'audioExamples/marimba/V2/loopbackFM.wav', ...
                 'audioExamples/marimba/V2/yLBFMMemb4.wav'};
    yMSWav = 'audioExamples/marimba/V2/modalSynthesis.wav';
    yMSPGWavs = {'audioExamples/marimba/V2/modalSynthesisPitchGlide.wav', ...
                'audioExamples/marimba/V2/yLBFMMemb4ModalSynthesisPitchGlide.wav'};
    resIRWav = '../proofOfConcept/resonatorIRs/marimbaTube.wav';
end

% Tom tom examples
if strcmp('tom tom', synthExample) == 1
    yLBFMWavs = {[baseDir 'audioExamples/tomtom/V2/loopbackFM.wav'], ...
                 [baseDir 'audioExamples/tomtom/V2/yLBFMMemb3.wav']};
    yMSWav = [baseDir 'audioExamples/tomtom/V2/modalSynthesis.wav'];
    yMSPGWavs = {[baseDir 'audioExamples/tomtom/V2/modalSynthesisPitchGlide.wav']...
                [baseDir 'audioExamples/tomtom/V2/yLBFMMemb3ModalSynthesisPitchGlide.wav']};
    resIRWav = '../proofOfConcept/resonatorIRs/taiko/taiko2cut.wav';
end

% Circular plate examples
if strcmp('circular plate', synthExample) == 1
    yLBFMWavs = {'audioExamples/circularPlate/V2/loopbackFM.wav', ...
                 'audioExamples/circularPlate/V2/yFBFMCP1.wav', ...
                 'audioExamples/circularPlate/V2/yFBFMCP2.wav', ...
                 'audioExamples/circularPlate/V2/ySAPFCP3.wav'};
    yMSWav = 'audioExamples/circularPlate/V2/modalSynthesis.wav';
    yMSPGWavs = {'audioExamples/circularPlate/V2/modalSynthesisPitchGlide.wav', ...
                 'audioExamples/circularPlate/V2/yFBFMCP1MSPG.wav', ...
                 'audioExamples/circularPlate/V2/yFBFMCP2MSPG.wav', ...
                 'audioExamples/circularPlate/V2/ySAPFCP3MSPG.wav'};
    resIRWav = 'resonatorIRs/CarpenterCenter.wav';
    %resIRWav = 'resonatorIRs/3000CStreetGarageStairwell.wav';
end

% Bass drum examples
% if strcmp('bass drum', synthExample) == 1
%     yLBFMWavs = {'../proofOfConcept/audioExamples/membrane/bass-drum/yFBFMMemb3.wav', ...
%                  '../proofOfConcept/audioExamples/membrane/bass-drum/ySAPFMemb3.wav'};
%     yMSWav = '../proofOfConcept/audioExamples/membrane/bass-drum/yMS.wav';
%     resIRWav = '../proofOfConcept/resonatorIRs/beduk/beduk2cut.wav';
% end

% Kettledrum examples
if strcmp('kettledrum', synthExample) == 1
    yLBFMWavs = {'audioExamples/kettledrum/V2/ySAPFMemb3.wav'};
    yMSWav = 'audioExamples/kettledrum/V2/modalSynthesis.wav';
    yMSPGWavs = {'audioExamples/kettledrum/V2/ySAPFMemb3MSPG.wav'};
    resIRWav = '../proofOfConcept/resonatorIRs/beduk/beduk2cut.wav';
end

% o-daiko examples
if strcmp('o-daiko', synthExample) == 1
    yLBFMWavs = {'audioExamples/odaiko/V2/yLBFMMemb3.wav'};
    yMSWav = 'audioExamples/odaiko/V2/modalSynthesis.wav';
    yMSPGWavs = {'audioExamples/odaiko/V2/yLBFMMemb3MSPG.wav'};
    resIRWav = '../proofOfConcept/resonatorIRs/beduk/beduk2cut.wav';
end

% Wood block
if strcmp('wood block', synthExample) == 1
    yLBFMWavs = {'audioExamples/woodBlock/V2/loopbackFM.wav'};
    yMSWav = 'audioExamples/woodBlock/V2/modalSynthesis.wav';
    yMSPGWavs = {'audioExamples/woodBlock/V2/modalSynthesisPitchGlide.wav'};
    resIRWav = '../proofOfConcept/resonatorIRs/113620__vidsyn__miscsoftnaturalgtrloud2-2.wav';
    %resIRWav = '../proofOfConcept/resonatorIRs/taiko/taiko2cut.wav';
end

[yMS, ~] = audioread(yMSWav);

addpath(genpath('../proofOfConcept'));


%% create excitations and perform commuted synthesis

excitations = zeros(N,2);
dexcitations = zeros(N,2);  % first derivative of excitations

%% generate raised cosine excitation

winLength = 8;

% raised cosine/Hann window
n = 0:winLength-1;
w = 0.5 * (1 - cos((2*pi*n)/(winLength-1)));

excitations(1:winLength, 1) = w;

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

audioOutDir = 'audioExamples/synthesisExampleSounds/V2/';
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
    
    
    if i==1
        % traditional modal synthesis with commuted synthesis
        rc_ms = percSynth(dexcitations(:,1), yMSWav, resIRWav);
        nb_ms = percSynth(dexcitations(:,2), yMSWav, resIRWav);
        audiowrite([audioOutDir synthExample '_MS_RC.wav'], scaleForSavingAudio(rc_ms(1:N)), fs);
        audiowrite([audioOutDir synthExample '_MS_NB.wav'], scaleForSavingAudio(nb_ms(1:N)), fs);
    end
    
    % write audio
    s = strsplit(yLBFMWavs{i}, '/');
    ss = strsplit(s{end}, '.');
    wavname = ss{1};
    
     % traditional modal synthesis with pitch glide with commuted synthesis
    rc_mspg = percSynth(dexcitations(:,1), yMSPGWavs{i}, resIRWav);
    nb_mspg = percSynth(dexcitations(:,2), yMSPGWavs{i}, resIRWav);
    audiowrite([audioOutDir synthExample '_' wavname '_LBFMMSPG_RC.wav'], scaleForSavingAudio(rc_mspg(1:N)), fs);
    audiowrite([audioOutDir synthExample '_' wavname '_LBFMMSPG_NB.wav'], scaleForSavingAudio(nb_mspg(1:N)), fs);
   
    
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

