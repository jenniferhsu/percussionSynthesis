% Chapter4Plots.m
% This script plots all the figures from Chapter 4.
%
% author: Jennifer Hsu
% date: Spring 2019

addpath(genpath('loopbackFMPercSynth/'));
savePlots = 1;
figDir = 'figures/';
filePrefix = 'Chapter4_';

%% Figure 4.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Musical Parameters / Static Timbre: 
% z_c,i oscillators or z_0,i oscillators, low f_c,i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 1;
f_low = 2000;

B = 0.9;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

N = fs*dur;

b0 = (sqrt(1-B^2) - 1)/B;

% === Modal synthesis frequencies ===
% generate modes for a 3x3 mode steel plate using the von Karman equations
modes = plateRectModes(3, 3, 1, 1);
modes = modes(:);
modes = unique(modes);

% f_low = 2000
fVecModal_low = modes*f_low;
fcVecModal_low = fVecModal_low./(sqrt(1-B^2));
Nf = length(fVecModal_low);

% === set up loopback FM MS structs/params ===
argStruct = struct();
argStruct.sinusoidArgs = struct();
argStruct.zcArgs = struct();
argStruct.z0Args = struct();

% set up envelope
e0Vec = ones(1, Nf)';
T60Vec = ones(1, Nf)';
env = envMat(e0Vec, T60Vec, dur, fs);

% === set up loopback FM MS arguments ===

argStruct.zcArgs.fcVec = fcVecModal_low;
argStruct.zcArgs.BVec = B*ones(1, Nf);
argStruct.zcArgs.BEndVec = B*ones(1, Nf);
argStruct.zcArgs.gVec = zeros(1, Nf);
argStruct.zcArgs.BGlideTypeVec = {'none', 'none', 'none'};

argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = fVecModal_low;
argStruct.z0Args.f0EndVec = fVecModal_low;
argStruct.z0Args.pitchGlideTypeVec = {'none', 'none', 'none'};
for f=1:Nf
    argStruct.z0Args.b0Mat(f,:) = b0*ones(1, N);
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f) = struct();
end

% === loopback FM synthesis ===
[yzcModal_low1, ~] = loopbackFMMS('zc', env, argStruct, fs);
[yz0Modal_low1, ~] = loopbackFMMS('z0', env, argStruct, fs);

% === plot ===
figure
subplot(211)
spectrogram(real(yzcModal_low1), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('Spectrogram of MS using z_{c,i}(n), f_{c,i}=2000Hz')
set(gca,'FontSize',15)
subplot(212)
spectrogram(real(yz0Modal_low1), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('Spectrogram of MS using z_{0,i}(n), f_{c,i}=2000Hz')
set(gca,'FontSize',15)
if savePlots
    saveas(gcf, [figDir filePrefix 'z0_zc_f_low'], 'epsc')
end

%% Figure 4.6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Musical Parameters / Static Timbre: 
% z_c,i oscillators or z_0,i oscillators, high f_c,i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 1;
f_high = 4000;

B = 0.9;    % feedback coefficient
g = 0.9999; % pitch glide coefficient
N = fs*dur;
b0 = (sqrt(1-B^2) - 1)/B;

% === Modal synthesis frequencies ===
% generate modes for a 3x3 mode steel plate using the von Karman equations
modes = plateRectModes(3, 3, 1, 1);
modes = modes(:);
modes = unique(modes);

% f_high = 5500;
fVecModal_high = modes*f_high;
fcVecModal_high = fVecModal_high./(sqrt(1-B^2));
Nf = length(fVecModal_low);

% === set up loopback FM MS structs/params ===
argStruct = struct();
argStruct.sinusoidArgs = struct();
argStruct.zcArgs = struct();
argStruct.z0Args = struct();

% set up envelope
e0Vec = ones(1, Nf)';
T60Vec = ones(1, Nf)';
env = envMat(e0Vec, T60Vec, dur, fs);

% === set up loopback FM MS arguments ===

argStruct.zcArgs.fcVec = fcVecModal_high;
argStruct.zcArgs.BVec = B*ones(1, Nf);
argStruct.zcArgs.BEndVec = B*ones(1, Nf);
argStruct.zcArgs.gVec = zeros(1, Nf);
argStruct.zcArgs.BGlideTypeVec = {'none', 'none', 'none'};

argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = fVecModal_high;
argStruct.z0Args.f0EndVec = fVecModal_high;
argStruct.z0Args.pitchGlideTypeVec = {'none', 'none', 'none'};
for f=1:Nf
    argStruct.z0Args.b0Mat(f,:) = b0*ones(1, N);
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f) = struct();
end

% === loopback FM synthesis ===
[yzcModal_high1, ~] = loopbackFMMS('zc', env, argStruct, fs);
[yz0Modal_high1, ~] = loopbackFMMS('z0', env, argStruct, fs);

% === plot ===
figure
subplot(211)
spectrogram(real(yzcModal_high1), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('Spectrogram of MS using z_{c,i}(n), f_{c,i}=4000Hz')
set(gca,'FontSize',15)
subplot(212)
spectrogram(real(yz0Modal_high1), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('Spectrogram of MS using z_{0,i}(n), f_{c,i}=4000Hz')
set(gca,'FontSize',15)
if savePlots
    saveas(gcf, [figDir filePrefix 'z0_zc_f_high'], 'epsc')
end

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

% === plot ===
figure
subplot(211)
spectrogram(real(mz0SameTimbre), hann(256), 128, 1024, fs, 'yaxis');
title('Loopback FM MS with b_{0,i}=0.6 for all modal components')
set(gca, 'FontSize', 15);
colorbar('off')
ylim([0 18])
subplot(212)
spectrogram(real(mz0DiffTimbre), hann(256), 128, 1024, fs, 'yaxis');
title('Loopback FM MS with b_{0,i}=[0.6, 0.3, 0.0] for the modal components')
set(gca, 'FontSize', 15);
colorbar('off')
ylim([0 18])
if savePlots
    saveas(gcf, [figDir filePrefix 'staticTimbreSameDiff'], 'epsc')
end

%% Figure 4.8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Musical Parameters / Time-varying Timbre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% === input parameters ===

% general
fs = 44100;
dur = 1;

% loopback FM parameters
f_low = 2000;
g = 0.9999; % pitch/timbre glide coefficient

% derived parameters
N = fs*dur;
BVec = g.^(0:N-1);             % feedback FM pitch glide coefficient
b0Vec = (sqrt(1-BVec.^2) - 1) ./ BVec;

% === Modal synthesis frequencies ===
% generate modes for a 3x3 mode steel plate using the von Karman equations

modes = plateRectModes(3, 3, 1, 1);
modes = modes(:);
modes = unique(modes);
Nf = length(modes);

% f_low = 2000
fVecModal_low = modes*f_low;
fcVecModal_low = fVecModal_low./(sqrt(1-BVec(end)^2));

% set up argument struct
argStruct = struct();
argStruct.sinusoidArgs = struct();
argStruct.zcArgs = struct();
argStruct.z0Args = struct();

% set up envelope
e0Vec = [1, 0.9, 0.89]';
T60Vec = [0.9, 0.88, 0.87]';
env = envMat(e0Vec, T60Vec, dur, fs);

% === set up loopback FM MS arguments ===

argStruct.zcArgs.fcVec = fcVecModal_low;
argStruct.zcArgs.BVec = B*ones(1, Nf);
argStruct.zcArgs.BEndVec = B*ones(1, Nf);
argStruct.zcArgs.gVec = zeros(1, Nf);
argStruct.zcArgs.BGlideTypeVec = {'useB', 'useB', 'useB'};
argStruct.zcArgs.BMat = repmat(BVec, Nf,1);

argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = fVecModal_low;
argStruct.z0Args.f0EndVec = fVecModal_low;
argStruct.z0Args.pitchGlideTypeVec = {'none', 'none', 'none'};
argStruct.z0Args.b0Mat = repmat(b0Vec, Nf, 1);
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f) = struct();
end

% === loopback FM synthesis ===

% time-varying timbre for z_c(n)
[zc_tv, ~] = loopbackFMMS('zc', env, argStruct, fs);

% time-varying timbre for z_0(n) 
[z0_tv, ~] = loopbackFMMS('z0', env, argStruct, fs);
stretchedAPFSynthesis(fVecModal_low, b0Vec, env, fs, [], 'none');

% static timbre z_0(n)
for f=1:Nf
    argStruct.z0Args.b0Mat(f,:) = b0Vec(1024)*ones(1, N);
end
[z0_static, ~] = loopbackFMMS('z0', env, argStruct, fs);


% === plot ===

figure
subplot(311)
spectrogram(real(z0_static), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('MS using z_{0,i}(n) oscillators, static timbre, b_{0,i}=-0.6312')
set(gca,'FontSize',14)
subplot(312)
spectrogram(real(z0_tv), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('MS using z_{0,i}(n) oscillators, time-varying timbre, b_{0,i}(n)')
set(gca,'FontSize',14)
subplot(313)
spectrogram(real(zc_tv), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('MS using z_{c,i}(n) oscillators, time-varying timbre, B_i(n)')
set(gca,'FontSize',14)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 8];
if savePlots
    saveas(gcf, [figDir filePrefix 'timeVaryingTimbre'], 'epsc')
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
title('Loopback FM MS with linear b_{0,i}')
set(gca, 'FontSize', 15);
colorbar('off')
ylim([0 18])
subplot(212)
spectrogram(real(mz0DiffTimbre), hann(256), 128, 1024, fs, 'yaxis');
title('Loopback FM MS with linear b_{0,i} and delayed higher frequencies')
set(gca, 'FontSize', 15);
colorbar('off')
ylim([0 18])
if savePlots
    saveas(gcf, [figDir filePrefix 'timeVaryingTimbreDelayedHigherFreqs'], 'epsc')
end

%% Figure 4.10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Musical Parameters / Sounding Frequency Diagram
% this figure compares how the spectrogram changes if we set the sounding
% frequencies to the modal frequencies vs. if we set the carrier frequencies
% to the modal frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 1;
B = 0.9;    % feedback coefficient

% === Modal synthesis frequencies ===
% generate modes for a 3x3 mode steel plate using the von Karman equations

modes = plateRectModes(3, 3, 1, 1);
modes = modes(:);
modes = unique(modes);
Nf = length(modes);

f_high2 = 5000;
fVecModal_high2 = modes*f_high2;
fcVecModal_high2 = fVecModal_high2./(sqrt(1-B^2));

% set up argument struct
argStruct = struct();
argStruct.sinusoidArgs = struct();
argStruct.zcArgs = struct();
argStruct.z0Args = struct();

% set up envelope
e0Vec = [1, 0.9, 0.89]';
T60Vec = [0.9, 0.88, 0.87]';
env = envMat(e0Vec, T60Vec, dur, fs);

% === set up loopback FM MS arguments ===

argStruct.zcArgs.fcVec = fVecModal_high2;
argStruct.zcArgs.BVec = B*ones(1, Nf);
argStruct.zcArgs.BEndVec = B*ones(1, Nf);
argStruct.zcArgs.gVec = zeros(1, Nf);
argStruct.zcArgs.BGlideTypeVec = {'none', 'none', 'none'};


% === loopback FM synthesis ===

% loopback FM with pitch glide with carrier frequencies = modal frequencies
[zc_fcIsModal, ~] = loopbackFMMS('zc', env, argStruct, fs);

% loopback FM with pitch glide with sounding frequencies = modal frequencies
% adjust the modal frequencies
argStruct.zcArgs.fcVec = fcVecModal_high2;
[zc_f0IsModal, ~] = loopbackFMMS('zc', env, argStruct, fs);


% === plots ===
figure
subplot(211)
spectrogram(real(zc_fcIsModal), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('Loopback FM modal synthesis with f_{c,i} = f_{r,i}')
set(gca,'FontSize',15)
subplot(212)
spectrogram(real(zc_f0IsModal), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('Loopback FM modal synthesis with f_{0,i} = f_{r,i}')
set(gca,'FontSize',15)
if savePlots
    saveas(gcf, [figDir filePrefix 'carrierOrSoundingAsModalFreq'], 'epsc')
end

%% Figure 4.11
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
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f).T60 = 1.0;
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
% Musical Parameters / Pitch Glide: 
% z_{c,i}(n) and z_{0,i}(n) create different spectrograms
% at high carrier frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 1;
f_high = 4000;

B = 0.9;    % feedback coefficient
g = 0.9999; % pitch glide coefficient
BVec = g.^(0:N-1);             % feedback FM pitch glide coefficient
N = fs*dur;
b0 = (sqrt(1-B^2) - 1)/B;
b0Vec = (sqrt(1-BVec.^2) - 1)./BVec;

% === Modal synthesis frequencies ===
% generate modes for a 3x3 mode steel plate using the von Karman equations
modes = plateRectModes(3, 3, 1, 1);
modes = modes(:);
modes = unique(modes);

% f_high = 5500;
fVecModal_high = modes*f_high;
fcVecModal_high = fVecModal_high./(sqrt(1-B^2));
Nf = length(fVecModal_low);

% === set up loopback FM MS structs/params ===
argStruct = struct();
argStruct.sinusoidArgs = struct();
argStruct.zcArgs = struct();
argStruct.z0Args = struct();

% set up envelope
e0Vec = ones(1, Nf)';
T60Vec = ones(1, Nf)';
env = envMat(e0Vec, T60Vec, dur, fs);

% === set up loopback FM MS arguments ===

argStruct.zcArgs.fcVec = fcVecModal_high;
argStruct.zcArgs.BVec = B*ones(1, Nf);
argStruct.zcArgs.BEndVec = B*ones(1, Nf);
argStruct.zcArgs.gVec = zeros(1, Nf);
argStruct.zcArgs.BGlideTypeVec = {'useB', 'useB', 'useB'};
argStruct.zcArgs.BMat = repmat(BVec, Nf, 1);

argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = fVecModal_high;
argStruct.z0Args.f0EndVec = fcVecModal_high;
argStruct.z0Args.pitchGlideTypeVec = {'expB', 'expB', 'expB'};
for f=1:Nf
    argStruct.z0Args.b0Mat(f,:) = b0Vec;
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f) = struct();
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f).g = g;
    argStruct.z0Args.zcArgsVec(f).wc = 2*pi*fcVecModal_high(f);
end

% === loopback FM synthesis ===
[yzcModal_high2, ~] = loopbackFMMS('zc', env, argStruct, fs);
[yz0Modal_high2, ~] = loopbackFMMS('z0', env, argStruct, fs);

% === plot ===
figure
subplot(211)
spectrogram(real(yzcModal_high2), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('Spectrogram of MS using z_{c,i}(n) oscillators with pitch glide, f_{c,i}=4000Hz')
set(gca,'FontSize',15)
subplot(212)
spectrogram(real(yz0Modal_high2), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off')
title('Spectrogram of MS using z_{0,i}(n) oscillators with pitch glide, f_{c,i}=4000Hz')
set(gca,'FontSize',15)
if savePlots
    saveas(gcf, [figDir filePrefix 'z0_zc_pitchGlide_f_high'], 'epsc')
end


%% Figure 4.13
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
