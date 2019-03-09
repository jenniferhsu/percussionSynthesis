% This script shows the difference between linear and nonlinear plates for
% Bilbao 2005's Sound Synthesis for Nonlinear Plates model.  

addpath(genpath('/Users/jenniferhsu/Documents/ucsd/fall2017/percussion_synthesis/code/matlab/nonlinearPlate'));

% set up parameter structure
params = {};
params.fs = 44100;
params.dur = 1.0;
params.outXFrac = 0.9;
params.outYFrac = 0.9;
params.Lx = 0.4;
params.Ly = 0.4;
params.h = 0.005;
params.E = 2*10^11;
params.rho = 7.86*10^3;
params.v = 0.3;
params.c = sqrt(173);
params.sigma = 1.38;
params.b1 = 0.002;

%% Example 1: Linear steel plate with high velocity input

params.NL = 0;
params.v0 = 70;

yLin = NLVonKarmanSteelPlateFDM(params);
yLin = yLin/max(abs(yLin));
yLin = 0.95*yLin;

audiowrite('audioExamples/linear.wav', yLin, fs);

%% Example 2: Nonlinear steel plate with medium velocity input

params.NL = 1;
params.v0 = 70;

yNL = NLVonKarmanSteelPlateFDM(params);
yNL = yNL/max(abs(yNL));
yNL = 0.95*yNL;

audiowrite('audioExamples/nonlinear.wav', yNL, fs);


%% plots
figure
subplot(211)
spectrogram(real(yLin), hann(256), 128, 1024, params.fs, 'yaxis');
title('linear steel plate')
subplot(212)
spectrogram(real(yNL), hann(256), 128, 1024, params.fs, 'yaxis');
title('nonlinear steel plate')


%% FFT peak picking

% FFT parameters
fs = params.fs;
dur = params.dur;

N = fs*dur;
Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

%% Linear Plate Peak Picking

% FFT analysis
YLin = fft(yLin, Nfft);
YLinPos = YLin(1:Nfft/2+1);
YLinPosdB = 20*log10(abs(YLinPos)/max(abs(YLinPos)));

% peak picking
[pksLin, locsLin] = findpeaks(YLinPosdB, 'MinPeakHeight', -60, 'MinPeakWidth', 2);

% plot
figure
plot(faxis, YLinPosdB)
hold on
for i=1:length(locsLin)
    plot(faxis(locsLin(i)), pksLin(i), 'r*');
end
ylim([-60 0])

% linear modal freqs
fLin = faxis(locsLin);

save('fLinMat', 'fLin');

%% Nonlinear Plate Peak Picking

% FFT analysis
YNL = fft(yNL, Nfft);
YNLPos = YNL(1:Nfft/2+1);
YNLPosdB = 20*log10(abs(YNLPos)/max(abs(YNLPos)));

% peak picking
[pksNL, locsNL] = findpeaks(YNLPosdB, 'MinPeakHeight', -60, 'MinPeakWidth', 2);

% plot
figure
plot(faxis, YNLPosdB)
hold on
for i=1:length(locsNL)
    plot(faxis(locsNL(i)), pksNL(i), 'r*');
end
plot(faxis, YLinPosdB, 'g--')
% for i=1:length(locsLin)
%     plot([faxis(locsLin(i)) faxis(locsLin(i))], [-60, 0], 'g--')
% end
ylim([-60 0])

% linear modal freqs
fNL = faxis(locsNL)