fs = 44100;
dur = 1;

T = 1/fs;
N = dur*fs;
nT = T*linspace(0, N-1, N);

% directories for saving
figDir = 'figures/';
audioDir = 'audioExamples/';

% ENVELOPE
env = ones(1, N);
fadeSec = 0.01;
fadeSamps = fadeSec*fs;
env(1:fadeSamps) = linspace(0, 1, fadeSamps);
env(end-fadeSamps+1:end) = linspace(1, 0, fadeSamps);

%% SINUSOID
f0 = 220;
f1 = 1000;
fd = ((f1-f0)/dur)*(nT.^2)/2 + f0*nT;

xs = sin(2*pi*fd);

xs = xs .* env;
xs = 0.95 * (xs/max(abs(xs)));

spectrogram(xs, hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
colormap('pink');
ylim([0 3])
set(gca, 'fontsize', 14);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2];
print([figDir 'whyLoopbackFM_sinusoid.eps'], '-depsc', '-r0')

audiowrite([audioDir 'whyLoopbackFM_sinusoid.wav'], xs, fs);

%% FM

fc = 220;
fm = 440;
MVec = linspace(0, 20, N);

xfm = cos(2*pi*fc*nT + MVec.*sin(2*pi*fm*nT));

xfm = xfm .* env;
xfm = 0.95 * (xfm/max(abs(xfm)));

spectrogram(xfm, hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
colormap('pink');
ylim([0 10])
set(gca, 'fontsize', 14);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2];
print([figDir 'whyLoopbackFM_FM.eps'], '-depsc', '-r0')

audiowrite([audioDir 'whyLoopbackFM_FM.wav'], xfm, fs);


spectrogram(xfm, hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
colormap('pink');
ylim([0 10])
set(gca, 'fontsize', 14);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2];
print([figDir 'whyLoopbackFM_FM.eps'], '-depsc', '-r0')

audiowrite([audioDir 'whyLoopbackFM_FM.wav'], xfm, fs);


% with a pitch glide
xfm = cos(2*pi*fd + 1*sin(2*pi*fm*nT));

xfm = xfm .* env;
xfm = 0.95 * (xfm/max(abs(xfm)));

spectrogram(xfm, hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
colormap('pink');
ylim([0 10])
set(gca, 'fontsize', 14);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2];
print([figDir 'whyLoopbackFM_FM.eps'], '-depsc', '-r0')

audiowrite([audioDir 'whyLoopbackFM_FM_pitchGlide.wav'], xfm, fs);


%% FEEDBACK FM

BVec = linspace(0, 1, N);

xfbfm = zeros(1, N);
xfbfm(1) = 1;

for n=2:N
    %xfbfm(n) = sin(2*pi*fm*n*T + BVec(n)*xfbfm(n-1));
    xfbfm(n) = sin(2*pi*fd(n) + BVec(n)*xfbfm(n-1));
end

xfbfm = xfbfm .* env;
xfbfm = 0.95 * (xfbfm/max(abs(xfbfm)));

spectrogram(xfbfm, hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
colormap('pink');
ylim([0 10])
set(gca, 'fontsize', 14);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2];
print([figDir 'whyLoopbackFM_feedbackFM.eps'], '-depsc', '-r0')

audiowrite([audioDir 'whyLoopbackFM_feedbackFM.wav'], xfbfm, fs);

%% FEEDBACK AM

BVec = linspace(0, 0.999, N);

xfbam = zeros(1, N);
xfbam(1) = 1;

for n=2:N
    %xfbam(n) = cos(2*pi*fm*n*T)*(1 + BVec(n)*xfbam(n-1));
    xfbam(n) = cos(2*pi*fd(n))*(1 + BVec(n)*xfbam(n-1));
end

xfbam = xfbam .* env;
xfbam = 0.95 * (xfbam/max(abs(xfbam)));

spectrogram(xfbam, hann(256), 128, 1024, fs, 'yaxis');
colorbar('off');
colormap('pink');
ylim([0 10])
set(gca, 'fontsize', 14);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2];
print([figDir 'whyLoopbackFM_feedbackAM.eps'], '-depsc', '-r0')

audiowrite([audioDir 'whyLoopbackFM_feedbackAM.wav'], xfbam, fs);