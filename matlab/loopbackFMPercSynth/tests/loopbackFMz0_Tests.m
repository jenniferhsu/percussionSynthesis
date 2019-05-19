% loopbackFMz0_Tests.m: tests for loopbackFMz0.m file
%loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, varargin);

addpath(genpath('../'));

%% static pitch and timbre

f0 = 304;
b0 = -0.623;
dur = 1;
fs = 44100;
[z01, pitchGlide1] = loopbackFMz0(f0, 0, 'none', b0, dur, fs);


% FFT analysis
N = dur*fs;
Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

Z01 = fft(z01, Nfft);
Z01Pos = Z01(1:Nfft/2+1);


%% time-varying pitch - linear

f0 = 304;
f0End = 625;
b0 = -0.623;
dur = 1;
fs = 44100;
[z02 , pitchGlide2] = loopbackFMz0(f0, f0End, 'lin', b0, dur, fs);

%% time-varying linear pitch and linear timbre

f0 = 304;
f0End = 625;
dur = 1;
fs = 44100;
b0 = linspace(0.001, 0.99, fs*dur);
[z03 , pitchGlide3] = loopbackFMz0(f0, f0End, 'lin', b0, dur, fs);

%% time-varying linear pitch and exponential timbre

f0 = 304;
f0End = 625;
dur = 1;
fs = 44100;
b0 = 0.9999.^(linspace(0, fs*dur-1, fs*dur));
[z04 , pitchGlide4] = loopbackFMz0(f0, f0End, 'lin', b0, dur, fs);


%% time-varying pitch - exponential, static timbre

f0 = 304;
f0End = 104;
dur = 1;
fs = 44100;
b0 = 0.5;
aStruct = struct();
aStruct.T60 = 0.8;
[z05 , pitchGlide5] = loopbackFMz0(f0, f0End, 'exp', b0, dur, fs, aStruct);


%% time-varying pitch - sqrt, static timbre

f0 = 304;
f0End = 104;
dur = 1;
fs = 44100;
b0 = 0.5;
[z06 , pitchGlide6] = loopbackFMz0(f0, f0End, 'sqrt', b0, dur, fs);

%% time-varying pitch - linear B

wcStruct = struct();
wcStruct.B = 0.999;
wcStruct.BEnd = 0.3;
wcStruct.wc = 2*pi*500;
b0 = 0.5;
dur = 1;
fs = 44100;
[z07 , pitchGlide7] = loopbackFMz0(0, 0, 'linB', b0, dur, fs, wcStruct);

%% time-varying pitch - exponential B

wcStruct = struct();
wcStruct.g = 0.9999;
wcStruct.wc = 2*pi*500;
b0 = 0.5;
dur = 1;
fs = 44100;
[z08 , pitchGlide8] = loopbackFMz0(0, 0, 'expB', b0, dur, fs, wcStruct);


%% plots

figure
plot(faxis, abs(Z01Pos))
title('loopbackFMz0: static pitch and timbre')

figure
spectrogram(real(z02), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide2/1000, 'r', 'linewidth', 2)
title('loopbackFMz0: linear pitch glide, static timbre')
ylim([0 5])

figure
spectrogram(real(z03), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide3/1000, 'r', 'linewidth', 2)
title('loopbackFMz0: linear pitch glide, linear timbre')
ylim([0 5])

figure
spectrogram(real(z04), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide4/1000, 'r', 'linewidth', 2)
title('loopbackFMz0: linear pitch glide, exponential timbre')
ylim([0 5])

figure
spectrogram(real(z05), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide5/1000, 'r', 'linewidth', 2)
title('loopbackFMz0: exponential pitch glide, static timbre')
ylim([0 5])

figure
spectrogram(real(z06), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide6/1000, 'r', 'linewidth', 2)
title('loopbackFMz0: sqrt pitch glide, static timbre')
ylim([0 5])

figure
spectrogram(real(z07), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide7/1000, 'r', 'linewidth', 2)
title('loopbackFMz0: linear B pitch glide, static timbre')
ylim([0 5])

figure
spectrogram(real(z08), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide8/1000, 'r', 'linewidth', 2)
title('loopbackFMz0: exponential B pitch glide, static timbre')
ylim([0 5])