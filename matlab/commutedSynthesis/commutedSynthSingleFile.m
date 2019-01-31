% commutedSynthSingleFile.m
% this script performs commuted synthesis for a single percussion sound
% file

addpath(genpath('../proofOfConcept'));

%% input parameters
fs = 44100;
N = 2*fs;

% modal synthesis signal location 
yMSWav = '../proofOfConcept/audioExamples/membrane/kettledrum-ySAPFMemb2.wav';
%yMSWav = 'audioExamples/timeVaryingAPF/yTVAPFModal2.wav';

% resonating body impulse response location
%resIRWav = '../proofOfConcept/resonatorIRs/195790__klankbeeld__cinematic-boom-130730-06.wav';
%resIRWav = '../proofOfConcept/resonatorIRs/10543__batchku__hit-low-14-002.aiff';
%resIRWav = '../proofOfConcept/resonatorIRs/39303__the-semen-incident__pvc-tub-bass.wav';
%resIRWav = '../proofOfConcept/resonatorIRs/45512__tim-kahn__whipper-1.flac';
%resIRWav = '../proofOfConcept/resonatorIRs/45595__tim-kahn__futon-bd-1.aiff';
%resIRWav = '../proofOfConcept/resonatorIRs/45596__tim-kahn__futon-bd-2.aiff';
%resIRWav = '../proofOfConcept/resonatorIRs/snareShell.wav';
resIRWav = '../proofOfConcept/resonatorIRs/taiko/taiko2.wav';
%resIRWav = '../proofOfConcept/resonatorIRs/beduk/beduk2.wav';

%resIRWav = 'resonatorIRs/CarpenterCenter.wav';
resIRWav = 'resonatorIRs/3000CStreetGarageStairwell.wav';

% input percussion signal for comparison
[yMS, ~] = audioread(yMSWav);
[yresIR, ~] = audioread(resIRWav);


%% create excitations and perform convolutional synthesis

excitations = zeros(N,3);
dexcitations = zeros(N,3);  % first derivative of excitations

%% 1. impulse

excitations(1,1) = 1;

% take derivative for velocity
dexcitations(:,1) = [diff(excitations(:,1)); 0]; 

y1 = percSynth(dexcitations(:,1), yMSWav, resIRWav);
y1 = y1(1:N);

%% 2. raised cosine

winLength = 64;

% raised cosine/Hann window
n = winLength/2:winLength-1;
w = 0.5 * (1 - cos((2*pi*n)/(winLength-1)));

excitations(1:winLength/2, 2) = ones(winLength/2, 1);
excitations(winLength/2+1:winLength, 2) = w;

% take derivative for velocity
dexcitations(:,2) = [diff(excitations(:,2)); 0]; 

y2 = percSynth(dexcitations(:,2), yMSWav, resIRWav);
y2 = y2(1:N);


%% 3. filtered noise burst

durNB = 0.15;
lowFreq = 200;
highFreq = 2000;

 % noise burst
sampNB = ceil(durNB*fs);
excitations(:,3) = [(2*rand(1, sampNB))-1 zeros(1, N-sampNB)]';

[B, A] = butter(5, [lowFreq/(fs/2) highFreq/(fs/2)], 'bandpass');
%freqz(B,A)
excitations(:,3) = filter(B, A, excitations(:,3));

% take derivative for velocity
dexcitations(:,3) = [diff(excitations(:,3)); 0]; 

y3 = percSynth(dexcitations(:,3), yMSWav, resIRWav);
y3 = y3(1:N);



