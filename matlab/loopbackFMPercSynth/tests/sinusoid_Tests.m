% loopbackFMz0_Tests.m: tests for loopbackFMz0.m file
%loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, varargin);

addpath(genpath('../'));

%% static pitch

f0 = 304;
dur = 1;
fs = 44100;
[s1, pitchGlide1] = sinusoid(f0, 0, 'none', dur, fs);


%% time-varying pitch - linear

f0 = 304;
f0End = 625;
dur = 1;
fs = 44100;
[s2 , pitchGlide2] = sinusoid(f0, f0End, 'lin', dur, fs);


%% time-varying pitch - exponential

f0 = 304;
f0End = 104;
dur = 1;
fs = 44100;
[s3 , pitchGlide3] = sinusoid(f0, f0End, 'exp', dur, fs);


%% time-varying pitch - sqrt

f0 = 304;
f0End = 104;
dur = 1;
fs = 44100;
[s4 , pitchGlide4] = sinusoid(f0, f0End, 'sqrt', dur, fs);

%% time-varying pitch - linear B

wcStruct = struct();
wcStruct.B = 0.999;
wcStruct.BEnd = 0.3;
wcStruct.wc = 2*pi*500;
dur = 1;
fs = 44100;
[s5 , pitchGlide5] = sinusoid(0, 0, 'linB', dur, fs, wcStruct);

%% time-varying pitch - exponential B

wcStruct = struct();
wcStruct.g = 0.9999;
wcStruct.wc = 2*pi*500;
dur = 1;
fs = 44100;
[s6 , pitchGlide6] = sinusoid(0, 0, 'expB', dur, fs, wcStruct);


%% plots

figure
spectrogram(real(s1), hann(256), 128, 1024, fs, 'yaxis');
hold on
%plot(1000*linspace(0, dur, N), pitchGlide1/1000, 'r', 'linewidth', 2)
title('sinusoid: no pitch glide')
ylim([0 5])

figure
spectrogram(real(s2), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide2/1000, 'r', 'linewidth', 2)
title('sinusoid: linear pitch glide')
ylim([0 5])

figure
spectrogram(real(s3), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide3/1000, 'r', 'linewidth', 2)
title('sinusoid: exponential pitch glide')
ylim([0 5])

figure
spectrogram(real(s4), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide4/1000, 'r', 'linewidth', 2)
title('sinusoid: sqrt pitch glide')
ylim([0 5])

figure
spectrogram(real(s5), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide5/1000, 'r', 'linewidth', 2)
title('sinusoid: linear B pitch glide')
ylim([0 5])

figure
spectrogram(real(s6), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(1000*linspace(0, dur, N), pitchGlide6/1000, 'r', 'linewidth', 2)
title('sinusoid: exponential B pitch glide')
ylim([0 5])