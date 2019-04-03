addpath(genpath('../proofOfConcept'));

%% input parameters
fs = 44100;
N = 2*fs;

synthExample = 'marimba';

if strcmp('marimba', synthExample) == 1
    yLBFMWav = '../proofOfConcept/audioExamples/marimba/yFBFMMemb3.wav';
    yMSWav = '../proofOfConcept/audioExamples/marimba/yMS.wav';
    resIRWav = '../proofOfConcept/resonatorIRs/marimbaTube.wav';
end


[yLBFM, ~] = audioread(yLBFMWav);
[yMS, ~] = audioread(yMSWav);
[resIR, ~] = audioread(resIRWav);

%% generate excitations
excitations = zeros(N,2);
dexcitations = zeros(N,2);  % first derivative of excitations

%% generate raised cosine excitations

winLength = 8;

% full cosine
wFull = 0.5 * (1 - cos((2*pi*(0:winLength-1))/(winLength-1)));
excitations(1:winLength,1) = wFull;

% half cosine
wHalf = 0.5 * (1 - cos((2*pi*(winLength/2:winLength-1))/(winLength-1)));
excitations(1:winLength/2, 2) = ones(winLength/2, 1);
excitations(winLength/2+1:winLength, 2) = wHalf;

% take derivative for velocity
dexcitations(:,1) = [diff(excitations(:,1)); 0]; 
dexcitations(:,2) = [diff(excitations(:,2)); 0]; 


%% synthesis

% basic modal synthesis
yMSFull = percSynth(dexcitations(:,1), yMSWav, resIRWav);
yMSFull = yMSFull(1:N);

yMSHalf = percSynth(dexcitations(:,2), yMSWav, resIRWav);
yMSHalf = yMSHalf(1:N);

% LBFM synthesis
yLBFMFull = percSynth(dexcitations(:,1), yLBFMWav, resIRWav);
yLBFMFull = yLBFMFull(1:N);

yLBFMHalf = percSynth(dexcitations(:,2), yLBFMWav, resIRWav);
yLBFMHalf = yLBFMHalf(1:N);




