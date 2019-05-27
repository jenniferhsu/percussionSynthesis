% plots and sound examples for the SMC presentation
% z0 and zc comparison


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

%% PART 1:
% we can change pitch and number of harmonics independently of one another


%% loopback FM (zc) sound example for increasing B

wc = 2*pi*440;
BInc = linspace(0, 0.999, N);

zcBInc = zeros(1, N);
zcBInc(1) = 1;
for i=2:N
    zcBInc(i) = exp(1j*wc*T*(1 + BInc(i)*real(zcBInc(1,i-1)))) * zcBInc(1,i-1);
end

zcBInc = zcBInc .* env;
w0Inc = wc * sqrt(1 - BInc.^2);


%% closed form (z0) for same pitch glide but keep number of harmonics the same

b0 = 0.6;

B = BInc(1);
BEnd = BInc(end);
k = (BEnd - B)/N;
l = B;
BVec = k*n + l;
w0 = wc * sqrt(1 - BVec.^2);
pitchGlide = w0/(2*pi);
Theta0 = (wc*T/(2*k)) .* ...
    ((k*n + l) .* sqrt(1 - (k*n + l).^2) + asin(k*n + l));  

z0BInc = (b0 + exp(1j*Theta0)) ./ (1 + b0*exp(1j*Theta0));
z0BInc = z0BInc .* env;

    
%% closed form (z0) for same harmonic number variation, but no pitch glide  
        
b0 = (sqrt(1 - BInc.^2) - 1) ./ BInc;

z0BIncTimbre = (b0 + exp(1j*wc*n*T)) ./ (1 + b0.*exp(1j*wc*n*T));
z0BIncTimbre = z0BIncTimbre .* env;

%% plot


figure
spectrogram(real(z0BInc), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1000, N), (w0Inc/2*pi)/10000, 'r', 'linewidth', 2);
%title('z_0(n) with pitch glide and no timbre variation');
ylim([0 7]);
colorbar('off');
set(gca, 'fontsize', 15);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4];
print([figDir 'loopbackFM_z0PitchGlide.eps'], '-depsc', '-r0')


figure
spectrogram(real(z0BIncTimbre), hann(256), 128, 1024, fs, 'yaxis');
%title('z_0(n) timbre variation and no pitch glide');
ylim([0 7]);
colorbar('off');
set(gca, 'fontsize', 15);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4];
print([figDir 'loopbackFM_z0HarmonicNumberVariation.eps'], '-depsc', '-r0')

%% save audio

audiowrite([audioDir 'loopbackFM_z0PitchGlide.wav'], scaleForSavingAudio(real(z0BInc)), fs);
audiowrite([audioDir 'loopbackFM_z0HarmonicNumberVariation.wav'], scaleForSavingAudio(real(z0BIncTimbre)), fs);


%% PART 2:
% they are the same, but loopback FM zc has some crazy aliasing due to
% numerical errors


%% loopback FM (zc) sound example for increasing B [wc high]

wc = 2*pi*8000;
B = 0.99;

zc_wcHigh = zeros(1, N);
zc_wcHigh(1) = 1;
for i=2:N
    zc_wcHigh(i) = exp(1j*wc*T*(1 + B*real(zc_wcHigh(1,i-1)))) * zc_wcHigh(1,i-1);
end

zc_wcHigh = zc_wcHigh .* env;


%% closed form (z0) for same pitch glide but keep timbre the same [wc high]

b0 = (sqrt(1 - B^2) - 1)/B;
w0 = wc * sqrt(1 - B^2);  

z0_wcHigh = (b0 + exp(1j*w0*n*T)) ./ (1 + b0*exp(1j*w0*n*T));
z0_wcHigh = z0_wcHigh .* env;




%% loopback FM (zc) sound example for increasing B [wc low]

wc = 2*pi*2000;
B = 0.99;

zc_wcLow = zeros(1, N);
zc_wcLow(1) = 1;
for i=2:N
    zc_wcLow(i) = exp(1j*wc*T*(1 + B*real(zc_wcLow(1,i-1)))) * zc_wcLow(1,i-1);
end

zc_wcLow = zc_wcLow .* env;


%% closed form (z0) for same pitch glide but keep timbre the same [wc low]

b0 = (sqrt(1 - B^2) - 1)/B;
w0 = wc * sqrt(1 - B^2);  

z0_wcLow = (b0 + exp(1j*w0*n*T)) ./ (1 + b0*exp(1j*w0*n*T));
z0_wcLow = z0_wcLow .* env;


%% loopback FM (zc) sound example for increasing B [wc high 2]

wc = 2*pi*12000;
B = 0.99;

zc_wcHigh2 = zeros(1, N);
zc_wcHigh2(1) = 1;
for i=2:N
    zc_wcHigh2(i) = exp(1j*wc*T*(1 + B*real(zc_wcHigh2(1,i-1)))) * zc_wcHigh2(1,i-1);
end

zc_wcHigh2 = zc_wcHigh2 .* env;


%% closed form (z0) for same pitch glide but keep timbre the same [wc high 2]

b0 = (sqrt(1 - B^2) - 1)/B;
w0 = wc * sqrt(1 - B^2);  

z0_wcHigh2 = (b0 + exp(1j*w0*n*T)) ./ (1 + b0*exp(1j*w0*n*T));
z0_wcHigh2 = z0_wcHigh2 .* env;

%% FFT analysis


Z0_wcLow = fft(z0_wcLow, Nfft);
Z0_wcHigh = fft(z0_wcHigh, Nfft);
Z0_wcHigh2 = fft(z0_wcHigh2, Nfft);

Zc_wcLow = fft(zc_wcLow, Nfft);
Zc_wcHigh = fft(zc_wcHigh, Nfft);
Zc_wcHigh2 = fft(zc_wcHigh2, Nfft);

Z0_wcLowPos = Z0_wcLow(1:Nfft/2+1);
Z0_wcHighPos = Z0_wcHigh(1:Nfft/2+1);
Z0_wcHigh2Pos = Z0_wcHigh2(1:Nfft/2+1);

Zc_wcLowPos = Zc_wcLow(1:Nfft/2+1);
Zc_wcHighPos = Zc_wcHigh(1:Nfft/2+1);
Zc_wcHigh2Pos = Zc_wcHigh2(1:Nfft/2+1);



%% plots


figure
%spectrogram(real(z0_wcHigh), hann(256), 128, 1024, fs, 'yaxis');
plot(faxis, 20*log10(abs(Z0_wcHighPos)/max(abs(Z0_wcHighPos))), 'linewidth', 2);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
%title('z_0(n), \omega _c = 8000');
xlim([0, 15000])
ylim([-60 0]);
set(gca, 'fontsize', 20);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
print([figDir 'loopbackFM_z0_wcHigh.eps'], '-depsc', '-r0')


figure
%spectrogram(real(zc_wcHigh), hann(256), 128, 1024, fs, 'yaxis');
plot(faxis, 20*log10(abs(Zc_wcHighPos)/max(abs(Zc_wcHighPos))), 'linewidth', 2);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
%title('z_c(n), \omega _c = 8000');
xlim([0, 15000])
ylim([-60 0]);
set(gca, 'fontsize', 20);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
print([figDir 'loopbackFM_zc_wcHigh.eps'], '-depsc', '-r0')

figure
%spectrogram(real(z0_wcLow), hann(256), 128, 1024, fs, 'yaxis');
plot(faxis, 20*log10(abs(Z0_wcLowPos)/max(abs(Z0_wcLowPos))), 'linewidth', 2);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
%title('z_0(n), \omega _c = 2000');
xlim([0, 15000])
ylim([-60 0]);
set(gca, 'fontsize', 20);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
print([figDir 'loopbackFM_z0_wcLow.eps'], '-depsc', '-r0')

figure
%spectrogram(real(zc_wcLow), hann(256), 128, 1024, fs, 'yaxis');
plot(faxis, 20*log10(abs(Zc_wcLowPos)/max(abs(Zc_wcLowPos))), 'linewidth', 2);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
%title('z_c(n), \omega _c = 2000');
xlim([0, 15000])
ylim([-60 0]);
set(gca, 'fontsize', 20);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
print([figDir 'loopbackFM_zc_wcLow.eps'], '-depsc', '-r0')

figure
%spectrogram(real(z0_wcHigh2), hann(256), 128, 1024, fs, 'yaxis');
plot(faxis, 20*log10(abs(Z0_wcHigh2Pos)/max(abs(Z0_wcHigh2Pos))), 'linewidth', 2);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
%title('z_c(n), \omega _c = 12000');
xlim([0, 15000])
ylim([-60 0]);
set(gca, 'fontsize', 20);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
print([figDir 'loopbackFM_z0_wcHigh2.eps'], '-depsc', '-r0')

figure
%spectrogram(real(zc_wcHigh2), hann(256), 128, 1024, fs, 'yaxis');
plot(faxis, 20*log10(abs(Zc_wcHigh2Pos)/max(abs(Zc_wcHigh2Pos))), 'linewidth', 2);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
%title('z_c(n), \omega _c = 12000');
xlim([0, 15000])
ylim([-60 0]);
set(gca, 'fontsize', 20);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
print([figDir 'loopbackFM_zc_wcHigh2.eps'], '-depsc', '-r0')

%% save audio examples

audiowrite([audioDir 'loopbackFM_z0_wcHigh2.wav'], scaleForSavingAudio(real(z0_wcHigh2)), fs);
audiowrite([audioDir 'loopbackFM_zc_wcHigh2.wav'], scaleForSavingAudio(real(zc_wcHigh2)), fs);
audiowrite([audioDir 'loopbackFM_z0_wcHigh.wav'], scaleForSavingAudio(real(z0_wcHigh)), fs);
audiowrite([audioDir 'loopbackFM_zc_wcHigh.wav'], scaleForSavingAudio(real(zc_wcHigh)), fs);
audiowrite([audioDir 'loopbackFM_z0_wcLow.wav'], scaleForSavingAudio(real(z0_wcLow)), fs);
audiowrite([audioDir 'loopbackFM_zc_wcLow.wav'], scaleForSavingAudio(real(zc_wcLow)), fs);