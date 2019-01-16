% percSynthTests.m
% this script shows how to use the percSynth function with different
% impulses

%% input parameters
fs = 44100;
N = 2*fs;
outputDir = 'audioExamples/convolutionalSynth/';

% modal synthesis signal and resonating body impulse response locations
yMSWav = 'audioExamples/yTVAPFModal2.wav';
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


%% 3. half a cycle of a short sinusoid
% decreasing sineLength retains more higher frequencies while increasing
% sineLength up to 4096 removes the higher frequencies
% sineLength=30 produces a cool industrial-type hit for TVAPF 2 files
%sineLength = 4096;
sineLength = 30;
halfsine = sin(2*pi*(1:sineLength/2)/sineLength); % half cycle of a sinusoid
excitations(1:sineLength/2, 3) = halfsine;
y3 = percSynth(excitations(:, 3), yMSWav, resIRwav);
y3 = y3(1:N);


%% 4. noise burst
durNB = .10;
sampNB = ceil(durNB*fs);
excitations(:,4) = [(2*rand(1, sampNB))-1 zeros(1, N-sampNB)];
y4 = percSynth(excitations(:,4), yMSWav, resIRwav);
y4 = y4(1:N);


%% 5. filtered noise burst
[B, A] = butter(5, [400/(fs/2) 2000/(fs/2)], 'bandpass');
%freqz(B,A)
excitations(:,5) = filter(B, A, excitations(:,4));
y5 = percSynth(excitations(:,5), yMSWav, resIRwav);
y5 = y5(1:N);


%% variable width pulse (hann window)
durPulse = .10; % this doesn't change the resulting sound really
sampPulse = ceil(durPulse*fs);
if mod(sampPulse, 2) ~= 0
    sampPulse = sampPulse + 1;
end
sampAttack = 50; % shortening this creates a brighter sound
w = hann(sampPulse)';
excitations(:,6) = [linspace(0, 1, sampAttack) ones(1, (sampPulse/2)-sampAttack) w(sampPulse/2:end) zeros(1, N-sampPulse-1)];
y6 = percSynth(excitations(:,6), yMSWav, resIRwav);
y6 = y6(1:N);


%% write to file
if ~exist(outputDir)
    mkdir(outputDir)
end

audiowrite([outputDir yMSWav(15:end-4) '_impulse.wav'], scaleForSavingAudio(y1), fs);
audiowrite([outputDir yMSWav(15:end-4) '_ASR.wav'], scaleForSavingAudio(y2), fs);
audiowrite([outputDir yMSWav(15:end-4) '_halfSine.wav'], scaleForSavingAudio(y3), fs);
audiowrite([outputDir yMSWav(15:end-4) '_noiseBurst.wav'], scaleForSavingAudio(y4), fs);
audiowrite([outputDir yMSWav(15:end-4) '_noiseBurstFiltered.wav'], scaleForSavingAudio(y5), fs);
audiowrite([outputDir yMSWav(15:end-4) '_varWidthPulse.wav'], scaleForSavingAudio(y6), fs);


%% plot the excitation functions used
figure
subplot(311)
plot((1:N)/fs, excitations(:,1), 'linewidth', 2);
xlabel('time (sec)');
ylabel('amplitude (linear)');
title('impulse');
grid on
subplot(312)
plot((1:N)/fs, excitations(:,2), 'linewidth', 2);
xlabel('time (sec)');
ylabel('amplitude (linear)');
title('adsr-type');
grid on
subplot(313)
plot((1:N)/fs, excitations(:,3), 'linewidth', 2);
xlabel('time (sec)');
ylabel('amplitude (linear)');
title('half sinusoid cycle');
grid on

%saveas(gcf, 'excitation_examples_pt1.png', 'png');

figure
subplot(311)
plot((1:N)/fs, excitations(:,4), 'linewidth', 2);
xlabel('time (sec)');
ylabel('amplitude (linear)');
title('noise burst');
grid on
subplot(312)
plot((1:N)/fs, excitations(:,5), 'linewidth', 2);
xlabel('time (sec)');
ylabel('amplitude (linear)');
title('filtered noise burst');
grid on
subplot(313)
plot((1:N)/fs, excitations(:,6), 'linewidth', 2);
xlabel('time (sec)');
ylabel('amplitude (linear)');
title('variable width pulse');
grid on

%saveas(gcf, 'excitation_examples_pt2.png', 'png');