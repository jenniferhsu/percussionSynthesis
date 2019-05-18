% Chapter4Plots.m
% This script plots all the figures from Chapter 4.

addpath(genpath('loopbackFMPercSynth/'));
savePlots = 1;
figDir = 'figures/';
filePrefix = 'Chapter4_';

%% Figure 4.7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Musical Parameters / Static Timbre: 
% all components have the same b0 or all components have different b0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 44100;
dur = 1.0;

% set up argument struct
argStruct = struct();
argStruct.sinusoidArgs = struct();
argStruct.zcArgs = struct();
argStruct.z0Args = struct();

% set up envelope
e0Vec = [1, 0.9, 0.89]';
T60Vec = [0.9, 0.88, 0.87]';
env = envMat(e0Vec, T60Vec, dur, fs);

Nf = length(e0Vec);
N = dur*fs;

% create the signal with the same b0
argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = [100, 1000, 4000];
argStruct.z0Args.f0EndVec = [100, 1000, 4000];
argStruct.z0Args.pitchGlideTypeVec = {'none', 'none', 'none'};
for f=1:Nf
    argStruct.z0Args.b0Mat(f,:) = 0.6*ones(1, N);
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f) = struct();
end
[mz0SameTimbre, ~] = loopbackFMMS('z0', env, argStruct, fs);

% create the signal with different types of pitch glides
for f=1:Nf
    argStruct.z0Args.b0Mat(f,:) = 0.9 * (Nf-f)/Nf*ones(1, N);
end
[mz0DiffTimbre, ~] = loopbackFMMS('z0', env, argStruct, fs);

% plot
figure
subplot(211)
spectrogram(real(mz0SameTimbre), hann(256), 128, 1024, fs, 'yaxis');
title('Loopback FM MS with b_0=0.6 for all modal components')
set(gca, 'FontSize', 15);
colorbar('off')
ylim([0 18])
subplot(212)
spectrogram(real(mz0DiffTimbre), hann(256), 128, 1024, fs, 'yaxis');
title('Loopback FM MS with b_0=[0.6, 0.3, 0.0] for the modal components')
set(gca, 'FontSize', 15);
colorbar('off')
ylim([0 18])
if savePlots
    saveas(gcf, [figDir filePrefix 'staticTimbreSameDiff'], 'epsc')
end

%% Figure 4.9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Musical Parameters / Time-varying Timbre: 
% alternate between linear and exponential timbre functions for b0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
T = 1/fs;
n = linspace(0, N-1, N);

% create the signal with linear timbre changes, all starting at the same
% time
argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = [663, 1616, 1987, 4780, 7749];
argStruct.z0Args.f0EndVec = [663, 1616, 1987, 4780, 7749];
argStruct.z0Args.pitchGlideTypeVec = {'none', 'none', 'none', 'none', 'none'};
for f=1:Nf
    argStruct.z0Args.b0Mat(f,:) = linspace(0.999*((Nf-(f-1))/Nf), 0.001, N);
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f) = struct();
end
[mz0SameTimbre, ~] = loopbackFMMS('z0', env, argStruct, fs);

% create the signal with linear timbre variation, but add a delay in both
% the timbre function and the envelope so that higher frequencies show up
% later
tau_d = -(N/2)*T/log(0.001/1);
d = 1 * exp(-n*T/tau_d);
for f=1:Nf
    NZeros = 1 + ((f-1)/Nf * 2000);
    NInc = 100;
    argStruct.z0Args.b0Mat(f,:) = [zeros(1, NZeros), ...
                                   linspace(0, 0.999, NInc), ...
                                   linspace(0.999*((Nf-(f-1))/Nf),0, N-NZeros-NInc)];
    env(f,1:NZeros) = zeros(1, NZeros);
    env(f,NZeros:(NZeros+NInc-1)) = linspace(0, e0Vec(f), NInc);
end
[mz0DiffTimbre, ~] = loopbackFMMS('z0', env, argStruct, fs);

% plot
figure
subplot(211)
spectrogram(real(mz0SameTimbre), hann(256), 128, 1024, fs, 'yaxis');
title('Loopback FM MS with linear b_0')
set(gca, 'FontSize', 15);
colorbar('off')
ylim([0 18])
subplot(212)
spectrogram(real(mz0DiffTimbre), hann(256), 128, 1024, fs, 'yaxis');
title('Loopback FM MS with linear b_0 and delayed higher frequencies')
set(gca, 'FontSize', 15);
colorbar('off')
ylim([0 18])
if savePlots
    saveas(gcf, [figDir filePrefix 'timeVaryingTimbreDelayedHigherFreqs'], 'epsc')
end


%% Figure 4.10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Musical Parameters / Time-Varying Pitch Glides: 
% all components have the same kind of pitch glide or the pitch glides for
% the different components have different kinds of pitch glides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% create the signal with the same type of pitch glides
argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = [663, 1616, 1987, 4780, 7749];
argStruct.z0Args.f0EndVec = [400, 800, 1293, 2355, 3458];
argStruct.z0Args.pitchGlideTypeVec = {'lin', 'lin', 'lin', 'lin', 'lin'};
for f=1:Nf
    argStruct.z0Args.b0Mat(f,:) = linspace(0.999*((Nf-(f-1))/Nf), 0.001, N);
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f) = struct();
end
[mz0SamePitchGlide, mz0Mat] = loopbackFMMS('z0', env, argStruct, fs);

% create the signal with different types of pitch glides
argStruct.z0Args.pitchGlideTypeVec = {'exp', 'lin', 'exp', 'lin', 'exp'};
[mz0DiffPitchGlide, mz0Mat] = loopbackFMMS('z0', env, argStruct, fs);

% plot
figure
subplot(211)
spectrogram(real(mz0SamePitchGlide), hann(256), 128, 1024, fs, 'yaxis');
title('Loopback FM MS with all linear pitch glides')
set(gca, 'FontSize', 15);
colorbar('off')
ylim([0 18])
subplot(212)
spectrogram(real(mz0DiffPitchGlide), hann(256), 128, 1024, fs, 'yaxis');
title('Loopback FM MS with alternating linear/exponential pitch glides')
set(gca, 'FontSize', 15);
colorbar('off')
ylim([0 18])
if savePlots
    saveas(gcf, [figDir filePrefix 'timeVaryingPitchGlideSameDiff'], 'epsc')
end


%% Figure 4.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Musical Parameters / Decay Time: 
% 5 modal frequencies that decay with different T60s and have different
% initial amplitudes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e0Vec = [1, 0.85, 0.7, 0.65 0.44]';
T60Vec = [0.98, 0.87, 0.73, 0.66, 0.45]';
fs = 44100;
dur = 1.0;
env = envMat(e0Vec, T60Vec, dur, fs);

Nf = size(env, 1);
N = size(env, 2);
nT = 0:1/fs:(dur-1/fs);

figure
for f=1:Nf
    plot(nT, env(f,:), 'linewidth', 2);
    hold on
end
xlabel('Time (seconds)');
ylabel('Amplitude (linear)');
title('Exponentially decaying amplitude envelopes')
grid on
set(gca, 'FontSize', 15);
if savePlots
    saveas(gcf, [figDir filePrefix 'decayTimeT60Envelopes'], 'epsc')
end
