% kick examples
addpath(genpath('../../loopbackFMPercSynth/'));
addpath(genpath('../../helperFunctions/'));
savePlots = 1;
saveAudio = 1;

figDir = 'figures/';
audioDir = 'audioExamples/';

%%
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

% loopback FM signal (z0)
fs = 44100;
dur = 1;
f0 = 100;
f0End = 40;
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

% save audio
if saveAudio
    audiowrite([audioDir 'hsu_kick_pitchGlidez0' '.wav'], ...
                scaleForSavingAudio(real(kick1)), fs);
end


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Exponentially decreasing pitch glide kick drum synthesis
% and linearly decreasing timbre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loopback FM signal (z0)
fs = 44100;
dur = 1;
f0 = 100;
f0End = 40;
N = fs*dur;
b0 = linspace(0.5, 0, N);
pitchGlideType = 'exp';
aStruct = struct();
aStruct.T60 = 0.6;
kick2 = loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, aStruct);

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
kick2 = kick2 .* w;


if saveAudio
    audiowrite([audioDir 'hsu_kick_pitchAndTimbreVariation' '.wav'], ...
                scaleForSavingAudio(real(kick2)), fs);
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Applying time-varying allpass filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loopback FM signal (z0)
fs = 44100;
dur = 1;
f0 = 100;
f0End = 40;
N = fs*dur;
b0 = linspace(0.5, 0, N);
pitchGlideType = 'exp';
aStruct = struct();
aStruct.T60 = 0.6;
kick3 = loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, aStruct);

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

% time-varying allpass filter parameters
TVAPFParams.fmVec = [100];
TVAPFParams.MVec = [120];
TVAPFParams.fbVec = [2500];
TVAPFParams.fpiVec = [f0];

[kick3, yMat] = applyTimeVaryingAPF2(kick3, w, fs, TVAPFParams);

if saveAudio
    audiowrite([audioDir 'hsu_kick_timeVaryingAPF' '.wav'], ...
                scaleForSavingAudio(real(kick3)), fs);
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kick Drum / Applying Commuted Synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resIRWav = '../../loopbackFMPercSynth/resonatorIRs/taiko/taiko2_fixedOffset.wav';
excitationType = 'rc';
excitationParams = struct();
excitationParams.winLength = 4;
kick4 = applyCommutedSynthesis(kick3, resIRWav, excitationType, excitationParams, fs);

if saveAudio
    audiowrite([audioDir 'hsu_kick_commutedSynthesis' '.wav'], ...
                scaleForSavingAudio(real(kick4)), fs);
end

%% plot

figure
subplot(221)
spectrogram(real(kick1), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Pitch glide')
xlim([0 1000])
ylim([0 1])
set(gca, 'FontSize', 12);
set(gca, 'XTick', [200 400 600 800 1000], 'XTickLabel', [200 400 600 800 1000])

% plot
subplot(222)
spectrogram(real(kick2), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Pitch glide and timbre variation');
xlim([0 1000])
ylim([0 1])
set(gca, 'FontSize', 12);
set(gca, 'XTick', [200 400 600 800 1000], 'XTickLabel', [200 400 600 800 1000])

subplot(223)
spectrogram(real(kick3), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Time-varying APF');
xlim([0 1000])
ylim([0 1])
set(gca, 'FontSize', 12);
set(gca, 'XTick', [200 400 600 800 1000], 'XTickLabel', [200 400 600 800 1000])

subplot(224)
spectrogram(real(kick4), hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
title('Time-varying APF and CS');
ylim([0 1])
set(gca, 'FontSize', 12);

sgtitle('Kick Drum');

if savePlots
    saveas(gcf, [figDir 'kickExample'], 'epsc')
end