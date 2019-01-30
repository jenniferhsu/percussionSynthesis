% percSynthTests.m
% this script shows how to use the percSynth function with different
% impulses

%% input parameters
fs = 44100;
N = 2*fs;
WRITETOFILE = 0; % set to 1 to save audio examples
outputDir = 'audioExamples/timeVaryingAPF/convolutionalSynth/';

% modal synthesis signal
%yMSWav = 'audioExamples/timeVaryingAPF/yTVAPFModal2.wav';
yMSWav = 'audioExamples/FBFMStretchedAPF/ySAPFModal1.wav';

% input percussion signal for comparison
[yMS, ~] = audioread(yMSWav);

%% create excitations and perform convolutional synthesis

excitations = zeros(N,6);
dexcitations = zeros(N,6); % first difference of excitations

%% 1. impulse
excitations(1,1) = 1;
dexcitations(:,1) = [diff(excitations(:,1)); 0];
y1 = percSynthExcitation(dexcitations(:,1), yMSWav);
y1 = y1(1:N);


%% 2. attack, sustain, release signal
% making attackSamps shorter retains more of the higher frequencies
% increasing attackSamps to 512 removes the higher frequencies and creates
% a very bassy sound.  attackSamps=8 sounds cool!

attackSamps = 2;
%attackSamps = 8;
%sustainSamps = 1024;
sustainSamps = 2;
releaseSamps = 1024;
excitations(1:attackSamps, 2) = linspace(0, 1, attackSamps);
excitations(attackSamps+1:attackSamps+1+sustainSamps, 2) = 1;
excitations(attackSamps+sustainSamps+1:attackSamps+sustainSamps+releaseSamps, 2) = linspace(1, 0, releaseSamps);

dexcitations(:,2) = [diff(excitations(:,2)); 0]; % take derivative for velocity?
y2 = percSynthExcitation(dexcitations(:,2), yMSWav);
y2 = y2(1:N);

% the part of the signal that looks and sounds good:
%soundsc(y2((attackSamps+sustainSamps+releaseSamps):end), fs)
y2 = y2/max(y2);
plot(y2((attackSamps+sustainSamps+releaseSamps):end))


%% 3. half a cycle of a short sinusoid
% decreasing sineLength retains more higher frequencies while increasing
% sineLength up to 4096 removes the higher frequencies
% sineLength=30 produces a cool industrial-type hit for TVAPF 2 files
%sineLength = 2048;
sineLength = 1024;
%sineLength = 30;
%sineLength = 4;
halfsine = sin(2*pi*(1:sineLength/2)/sineLength); % half cycle of a sinusoid
excitations(1:sineLength/2, 3) = halfsine;

dexcitations(:,3) = [diff(excitations(:,3)); 0]; % take derivative for velocity?
y3 = percSynthExcitation(dexcitations(:,3), yMSWav);

y3 = y3(1:N);

% the part of the signal that looks and sounds good:
%soundsc(y3(sineLength/2:end), fs)

%% 4. noise burst
durNB = .10;
%durNB = .2;
sampNB = ceil(durNB*fs);
excitations(:,4) = [(2*rand(1, sampNB))-1 zeros(1, N-sampNB)];

dexcitations(:,4) = [diff(excitations(:,4)); 0];
y4 = percSynthExcitation(dexcitations(:,4), yMSWav);
y4 = y4(1:N);


%% 5. filtered noise burst
[B, A] = butter(5, [400/(fs/2) 2000/(fs/2)], 'bandpass');
%freqz(B,A)
excitations(:,5) = filter(B, A, excitations(:,4));

dexcitations(:,5) = [diff(excitations(:,5)); 0];
y5 = percSynthExcitation(dexcitations(:,5), yMSWav);
y5 = y5(1:N);


%% variable width pulse (hann window)
durPulse = .10; 
durPulse = .001;

sampPulse = ceil(durPulse*fs);
if mod(sampPulse, 2) ~= 0
    sampPulse = sampPulse + 1;
end
%sampAttack = 50; % shortening this creates a brighter sound
w = hann(sampPulse)';
%excitations(:,6) = [linspace(0, 1, sampAttack) ones(1, (sampPulse/2)-sampAttack) w(sampPulse/2:end) zeros(1, N-sampPulse-1)];
excitations(:,6) = [w(sampPulse/2:end) zeros(1, N-sampPulse/2-1)];

dexcitations(:,6) = [diff(excitations(:,6)); 0];
y6 = percSynthExcitation(dexcitations(:,6), yMSWav);
y6 = y6(1:N);

% the part of the signal that looks and sounds good:
%soundsc(y3(sampPulse/2:end), fs)


%% write to file

if WRITETOFILE==1
    if ~exist(outputDir)
        mkdir(outputDir)
    end

    audiowrite([outputDir yMSWav(15:end-4) '_impulse.wav'], scaleForSavingAudio(y1), fs);
    audiowrite([outputDir yMSWav(15:end-4) '_ASR.wav'], scaleForSavingAudio(y2), fs);
    audiowrite([outputDir yMSWav(15:end-4) '_halfSine.wav'], scaleForSavingAudio(y3), fs);
    audiowrite([outputDir yMSWav(15:end-4) '_noiseBurst.wav'], scaleForSavingAudio(y4), fs);
    audiowrite([outputDir yMSWav(15:end-4) '_noiseBurstFiltered.wav'], scaleForSavingAudio(y5), fs);
    audiowrite([outputDir yMSWav(15:end-4) '_varWidthPulse.wav'], scaleForSavingAudio(y6), fs);
end

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