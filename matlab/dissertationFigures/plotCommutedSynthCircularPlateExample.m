% plotCommutedSynthCircularPlateExample.m
%
% plots the spectrogram of the circular plate synthesis with and without
% commuted synthesis

plotOutDir = 'figures/'

% read in files
CPnoCSWav = 'audioExamples/circularPlate/loopbackFM.wav';
CPCSWav = 'audioExamples/circularPlate/commutedSynth/loopbackFM.wav';

[CPnoCS, fs] = audioread(CPnoCSWav);
[CPCS, ~] = audioread(CPCSWav);


% plot spectrograms
figure
spectrogram(real(CPnoCS), hann(256), 128, 1024, fs, 'yaxis');
ylim([0 16])
colorbar('off')
fig = gcf
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 1.6];
print([plotOutDir 'commutedSynthExamplesNoCS'], '-depsc', '-r0')

figure
spectrogram(real(CPCS), hann(256), 128, 1024, fs, 'yaxis');
ylim([0 16])
colorbar('off')
fig = gcf
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 1.6];
print([plotOutDir 'commutedSynthExamplesWithCS'], '-depsc', '-r0')
