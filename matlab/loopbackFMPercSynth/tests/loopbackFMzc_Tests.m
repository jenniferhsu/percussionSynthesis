% loopbackFMzc_Tests.m: tests for loopbackFMzc.m file

addpath(genpath('../'));

%% static pitch and timbre

fc = 698;
B = 0.9;
dur = 1;
fs = 44100;
[zc1, pitchGlide1] = loopbackFMzc(fc, B, 0, 0, 'none', dur, fs);

% FFT analysis
N = dur*fs;
Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

ZC1 = fft(zc1, Nfft);
ZC1Pos = ZC1(1:Nfft/2+1);

% check w0
wc = 2*pi*fc;
w0 = wc * sqrt(1 - B^2);
f0 = w0/(2*pi)

%% time-varying pitch and timbre - use B

fc = 698;
g = 0.9999;
dur = 1;
fs = 44100;
N = fs*dur;
B = [linspace(0, 0.5, N/2) linspace(0.5, 0.2, N/2)];
[zc2 , pitchGlide2] = loopbackFMzc(fc, B, 0, g, 'useB', dur, fs);


%% time-varying pitch and timbre - exponential B

fc = 698;
g = 0.9999;
dur = 1;
fs = 44100;
[zc3 , pitchGlide3] = loopbackFMzc(fc, 0, 0, g, 'expB', dur, fs);

%% time-varying pitch and timbre - linear B

fc = 698;
B = 0.999;
BEnd = 0.01;
dur = 1;
fs = 44100;
[zc4, pitchGlide4] = loopbackFMzc(fc, B, BEnd, 0, 'linB', dur, fs);


%% plots

figure
plot(faxis, abs(ZC1Pos))
title('loopbackFMzc: static pitch and timbre')

figure
spectrogram(real(zc2), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide2/1000, 'r', 'linewidth', 2)
title('loopbackFMzc: useB pitch glide')
ylim([0 5])

figure
spectrogram(real(zc3), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide3/1000, 'r', 'linewidth', 2)
title('loopbackFMzc: expB pitch glide')
ylim([0 5])

figure
spectrogram(real(zc4), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide4/1000, 'r', 'linewidth', 2)
title('loopbackFMzc: linB pitch glide')
ylim([0 5])