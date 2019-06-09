% Chapter5Plots.m
% This script plots all the figures from Chapter 5.
%
% author: Jennifer Hsu
% date: Spring 2019

addpath(genpath('../loopbackFMPercSynth/'));
addpath(genpath('../helperFunctions/'));
savePlots = 0;
saveAudio = 1;

audioDir = 'audioExamples/';
filePrefix = 'Chapter5_';

%% Figure 5.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Static Pitch and Timbre: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loopback FM signal (z0)
fs = 44100;
dur = 1.0;
f0 = 43.654;
f0End = 100;
b0 = 0.2;
pitchGlideType = 'none';
kick0 = loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs);

% envelope
A0 = 1;
T60 = 0.8;

T = 1/fs;
N = dur*fs;
n = 0:N-1;

n60 = T60*fs;
A60 = A0/10^(60/20);
tau = -n60*T/log(A60/A0);

w = A0*exp(-n*T/tau);
kick0 = kick0 .* w;

% loopback FM (zc)
B = -2*b0/(b0^2 + 1);
BEnd = B;
BGlideType = 'none';
g = 0;
fc = f0/sqrt(1 - B^2);
kick0c = loopbackFMzc(fc, B, BEnd, g, BGlideType, dur, fs);
kick0c = kick0c .* w;

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_kick_staticPitchAndTimbre' '.wav'], ...
                scaleForSavingAudio(real(kick0)), fs);
end

%% Figure 5.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Exponentially decreasing pitch glide: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
dur = 1.0;
T = 1/fs;
N = dur*fs;
n = 0:N-1;

% exponentially decreasing pitch glide parameters
fx = 100;
fy = 40;
r = 0.001;
A = 1;
t_d = 0.6;

n_d = t_d * fs;

% calculate the exponentially decreasing pitch glide from 0 to 1
tau = -(n_d*T)/log(r/A);
d = exp(-n*T/tau);

% scale using fx and fy
f0n = (fx - fy) * d + fy;


%% Figure 5.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Exponentially decreasing pitch glide B(n) for z_{c,i}(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0n = 2*pi*f0n;
wc = w0n(1);
Bn = sqrt(1 - (w0n./wc).^2);


%% Figure 5.4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Exponentially decreasing pitch glide kick drum synthesis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loopback FM signal (z0)
fs = 44100;
dur = 1;
f0 = 100;
f0End = 43.654;
b0 = 0.2;
pitchGlideType = 'exp';
aStruct = struct();
aStruct.T60 = 0.6;
kick1 = loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, aStruct);

% envelope
A0 = 1;
T60 = 0.8;

T = 1/fs;
N = dur*fs;
n = 0:N-1;

n60 = T60*fs;
A60 = A0/10^(60/20);
tau = -n60*T/log(A60/A0);

w = A0*exp(-n*T/tau);
kick1 = kick1 .* w;

% loopback FM (zc)
w0n = 2*pi*f0n;
wc = w0n(1);
Bn = sqrt(1 - (w0n./wc).^2);

kick1c = loopbackFMzc(fc, Bn, 0, 0, 'useB', dur, fs);
kick1c = kick1c .* w;

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_kick_pitchGlidez0' '.wav'], ...
                scaleForSavingAudio(real(kick1)), fs);
    audiowrite([audioDir 'hsu_kick_pitchGlidezc' '.wav'], ...
                scaleForSavingAudio(real(kick1c)), fs);
end
