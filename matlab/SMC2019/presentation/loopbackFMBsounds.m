% plots and sound examples for the SMC presentation
% how does changing B affect the sample-by-sample rotation equation
% results?

addpath(genpath('../../proofOfConcept'))

fs = 44100;
dur = 1.0;
T = 1/fs;

N = fs*dur;
n = linspace(0, N-1, N);

% directories for saving
figDir = 'figures/';
audioDir = 'audioExamples/';

% fft parameters
Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

%% envelope
fadeTime = 0.1; % seconds
fadeSamps = fadeTime*fs;
env = [ones(1, N-fadeSamps) linspace(1, 0, fadeSamps)];


%% loopback FM sound examples for different B

wc = 2*pi*440;
BVals = [0:0.2:0.8 0.9 0.99];
NB = length(BVals);

w0Vals = wc*sqrt(1 - BVals.^2);
f0Vals = w0Vals/(2*pi);

zcMat = zeros(NB, N);
zcMat(:,1) = 1;

for nb=1:NB
    B = BVals(nb);
    for i=2:N
        zcMat(nb,i) = exp(1j*wc*T*(1 + B*real(zcMat(nb,i-1)))) * zcMat(nb,i-1);
    end
    zcMat(nb,:) = zcMat(nb,:) .* env;
end

%% loopback FM sound example for increasing B

BInc = linspace(0, 0.999, N);

zcBInc = zeros(1, N);
zcBInc(1) = 1;
for i=2:N
    zcBInc(i) = exp(1j*wc*T*(1 + BInc(i)*real(zcBInc(1,i-1)))) * zcBInc(1,i-1);
end

zcBInc = zcBInc .* env;
w0Inc = wc * sqrt(1 - BInc.^2);

%% loopback FM sound example for decreasing B

BDec = linspace(0.999, 0, N);

zcBDec = zeros(1, N);
zcBDec(1) = 1;
for i=2:N
    zcBDec(i) = exp(1j*wc*T*(1 + BDec(i)*real(zcBDec(1,i-1)))) * zcBDec(1,i-1);
end

zcBDec = zcBDec .* env;
w0Dec = wc * sqrt(1 - BDec.^2);

%% FFT analysis

ZcPosMat = zeros(NB, Nfft/2+1);
for nb=1:NB
    Zc = fft(zcMat(nb,:), Nfft);
    ZcPosMat(nb,:) = Zc(1:Nfft/2+1);
end


%% plots

for nb=1:NB
    B = BVals(nb);
    
    figure
    %spectrogram(real(zcMat(nb,:)), hann(256), 128, 1024, fs, 'yaxis');
    plot(faxis, 20*log10(abs(ZcPosMat(nb,:))/max(abs(ZcPosMat(nb,:)))), 'linewidth', 2);
    hold on
    plot([wc/(2*pi) wc/(2*pi)], [-60 0], 'g', 'linewidth', 2);
    plot([f0Vals(nb) f0Vals(nb)], [-60 0], 'r--', 'linewidth', 2);
    %title(sprintf('B = %.2f, \\omega _0 = 2\\pi%.1f ', BVals(nb), w0Vals(nb)/(2*pi)));
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    xlim([0 3000]);
    ylim([-60 0]);
    grid on
    set(gca, 'fontsize', 15);
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2];
    print([figDir 'loopbackFM_BVal' num2str(B) '.eps'], '-depsc', '-r0')
    close
end

% B increasing
figure
spectrogram(real(zcBInc), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1000, N), (w0Inc/2*pi)/10000, 'r', 'linewidth', 2);
%title('B increasing from 0 to 0.99');
ylim([0 3]);
colorbar('off');
set(gca, 'fontsize', 15);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 5];
print([figDir 'loopbackFM_BInc.eps'], '-depsc', '-r0')

% B decreasing
figure
spectrogram(real(zcBDec), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1000, N), (w0Dec/2*pi)/10000, 'r', 'linewidth', 2);
%title('B decreasing from 0.99 to 0.0');
ylim([0 3]);
colorbar('off');
set(gca, 'fontsize', 15);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 5];
print([figDir 'loopbackFM_BDec.eps'], '-depsc', '-r0')


%% B increasing plot #2 for z0 and zc comparison slide

figure
spectrogram(real(zcBInc), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1000, N), (w0Inc/2*pi)/10000, 'r', 'linewidth', 2);
%title('B increasing from 0 to 0.99');
ylim([0 7]);
colorbar('off');
set(gca, 'fontsize', 15);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4];
print([figDir 'loopbackFM_BInc_smaller.eps'], '-depsc', '-r0')

%% save audio

for nb=1:NB
    B = BVals(nb);
    audiowrite([audioDir 'loopbackFM_BVal' num2str(B) '.wav'], scaleForSavingAudio(real(zcMat(nb,:))), fs);
end

audiowrite([audioDir 'loopbackFM_BInc.wav'], scaleForSavingAudio(real(zcBInc)), fs);
audiowrite([audioDir 'loopbackFM_BDec.wav'], scaleForSavingAudio(real(zcBDec)), fs);

