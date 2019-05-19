% loopbackFMMS_Tests.m tests for loopbackFMMS.m file
%
% author: Jennifer Hsu
% date: Spring 2019

addpath(genpath('../'));
fs = 44100;
dur = 1.0;

% set up argument struct
argStruct = struct();
argStruct.sinusoidArgs = struct();
argStruct.zcArgs = struct();
argStruct.z0Args = struct();

% set up envelope
e0Vec = [1, 0.9, 0.89, 0.8, 0.75]';
T60Vec = [0.9, 0.88, 0.87, 0.85, 0.77]';
env = envMat(e0Vec, T60Vec, dur, fs);

Nf = length(e0Vec);
N = dur*fs;


%% traditional MS

argStruct.sinusoidArgs.f0Vec = [100, 2500, 3330, 5600, 6980];
argStruct.sinusoidArgs.f0EndVec = [70, 1250, 2110, 4300, 5550];
argStruct.sinusoidArgs.pitchGlideTypeVec = {'lin', 'exp', 'lin', 'exp', 'lin'};
argStruct.sinusoidArgs.zcArgsVec = [0, 0, 0, 0, 0];

[ms, msMat] = loopbackFMMS('s', env, argStruct, fs);


%% loopback FM zc MS

argStruct.zcArgs.fcVec = [4700, 5772, 6366, 7967, 8948];
argStruct.zcArgs.BVec = [0.99 0.96 0.95 0.8 0.5];
argStruct.zcArgs.BEndVec = [0.5 0.4 0.3 0.2 0.0];
argStruct.zcArgs.gVec = [0.9999, 0.999, 0, 0.99, 0];
argStruct.zcArgs.BGlideType = {'expB', 'expB', 'linB', 'expB', 'linB'};

[mzc, mzcMat] = loopbackFMMS('zc', env, argStruct, fs);

%% loopback FM z0 MS
argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = [663, 1616, 1987, 4780, 7749];
argStruct.z0Args.f0EndVec = [400, 800, 1293, 2355, 3458];
argStruct.z0Args.pitchGlideTypeVec = {'exp', 'lin', 'sqrt', 'linB', 'expB'};
for f=1:Nf
    argStruct.z0Args.b0Mat(f,:) = linspace(0.999*((Nf-(f-1))/Nf), 0.001, N);
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f) = struct();
end
argStruct.z0Args.zcArgsVec(4).B = 0.999;
argStruct.z0Args.zcArgsVec(4).BEnd = 0.5;
argStruct.z0Args.zcArgsVec(4).wc = 2*pi*2000;
argStruct.z0Args.zcArgsVec(5).g = 0.9999;
argStruct.z0Args.zcArgsVec(5).wc = 2*pi*1000;

[mz0, mz0Mat] = loopbackFMMS('z0', env, argStruct, fs);

%% plots

figure
spectrogram(real(ms), hann(256), 128, 1024, fs, 'yaxis');
title('traditional MS')
ylim([0 18])

figure
spectrogram(real(mzc), hann(256), 128, 1024, fs, 'yaxis');
title('loopback FM (zc) MS')
ylim([0 18])

figure
spectrogram(real(mz0), hann(256), 128, 1024, fs, 'yaxis');
title('loopback FM (z0) MS')
ylim([0 18])