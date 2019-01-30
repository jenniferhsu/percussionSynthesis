% percSynthTestNL.m
% this script compares percSynth.m and percSynthNL.m to see if there is a
% difference betwee the output

%% input parameters
fs = 44100;
N = 2*fs;
%outputDir = 'audioExamples/convolutionalSynth/';

% modal synthesis signal 
%yMSWav = 'audioExamples/timeVaryingAPF/yTVAPFModal2.wav';
%yMSWav = 'audioExamples/FBFMRotational/yFBFMMesh4.wav';
yMSWav = 'audioExamples/FBFMStretchedAPF/ySAPFMesh4.wav';

% resonating body impulse response locations
resIRwav = 'resonatorIRs/CarpenterCenter.wav';
%resIRwav = 'resonatorIRs/3000CStreetGarageStairwell.wav';

% input percussion signal for comparison
[yMS, ~] = audioread(yMSWav);


%% create excitations and perform convolutional synthesis

excitations = zeros(N,6);

%% 1. impulse
excitations(1,1) = 1;
y1 = percSynth(excitations(:,1), yMSWav, resIRwav);
y1 = y1(1:N);

y1NL = percSynthNL(excitations(:,1), yMSWav, resIRwav);
y1NL = y1NL(1:N);

sum(y1-y1NL)

figure
subplot(211)
plot(y1)
subplot(212)
plot(y1NL)

figure
subplot(211)
spectrogram(y1, hann(256), 128, 1024, fs, 'yaxis')
subplot(212)
spectrogram(y1NL, hann(256), 128, 1024, fs, 'yaxis')


%% 2. attack, sustain, release signal
% making attackSamps shorter retains more of the higher frequencies
% increasing attackSamps to 512 removes the higher frequencies and creates
% a very bassy sound.  attackSamps=8 sounds cool!
%attackSamps = 512;
attackSamps = 8;
sustainSamps = 1024;
releaseSamps = 8*attackSamps;
excitations(1:attackSamps, 2) = linspace(0, 1, attackSamps);
excitations(attackSamps+1:attackSamps+1+sustainSamps, 2) = 1;
excitations(attackSamps+sustainSamps+1:attackSamps+sustainSamps+releaseSamps, 2) = linspace(1, 0, releaseSamps);
y2 = percSynth(excitations(:,2), yMSWav, resIRwav);
y2 = y2(1:N);

y2NL = percSynth(excitations(:,2), yMSWav, resIRwav);
y2NL = y2NL(1:N);

sum(y2-y2NL)

figure
subplot(211)
plot(y2)
subplot(212)
plot(y2NL)

figure
subplot(211)
spectrogram(y2, hann(256), 128, 1024, fs, 'yaxis')
subplot(212)
spectrogram(y2NL, hann(256), 128, 1024, fs, 'yaxis')

