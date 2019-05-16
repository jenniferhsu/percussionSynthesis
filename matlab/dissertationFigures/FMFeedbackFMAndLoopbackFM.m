fs = 44100;
dur = 1;
saveDissertationFigures = 1;

% carrier and modulator frequencies
fc = 200;
fm = 150;
I = 1;
B = 0.9;

% derived parameters
N = fs*dur;
T = 1/fs;
nT = (0:N-1) * T;

% fft parameters
Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);


%% FM
%%%%%%%%%%%%%%%%%%%
zfmc = zeros(1, N);
zfmm = zeros(1, N);
zfmc(1) = 1;
zfmm(1) = 1;
%zfmm = I * exp(j*2*pi*fm*(0:N-1)*T);
for n=2:N
    % modulator
    zfmm(n) = exp(1j*2*pi*fm*T) * zfmm(n-1);
    % FM synthesis
    zfmc(n) = exp(1j*(2*pi*fc*T + imag(zfmm(n)-zfmm(n-1)))) * zfmc(n-1);
end

Zfmm = fft(zfmm, Nfft);
ZfmmPos = Zfmm(1:Nfft/2+1);

Zfmc = fft(zfmc, Nfft);
ZfmcPos = Zfmc(1:Nfft/2+1);

figure
plot(faxis, 20*log10(abs(ZfmmPos)/max(abs(ZfmmPos))));
hold on
plot([fm fm], [-60 0], 'r--')
ylim([-60 0]);
title('zfmm magnitude spectrum');

figure
plot(faxis, 20*log10(abs(ZfmcPos)/max(abs(ZfmcPos))));
hold on
plot([fc fc], [-60 0], 'r--')
% for i=-1:4
%     plot([fc+(i*fm) fc+(i*fm)], [-60 0], 'r--')
% end
ylim([-60 0]);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Magnitude spectrum of z_{fm,c}(n)');
set(gca, 'FontSize', 15)
grid on;
legend('z_{fm,c}', 'f_c');
if saveDissertationFigures==1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, 'figures/FMMagSpec', 'epsc')
end


xfm = cos(2*pi*fc*(0:N-1)*T + I*sin(2*pi*fm*(0:N-1)*T));
Xfm = fft(xfm, Nfft);
XfmPos = Xfm(1:Nfft/2+1);
figure
plot(faxis, 20*log10(abs(XfmPos)/max(abs(XfmPos))));
hold on
plot([fc fc], [-60 0], 'r--')
ylim([-60 1])


%% Feedback FM
%%%%%%%%%%%%%%%%%%%

zfb = zeros(1, N);
zfb(1) = 1;
zfb(2) = 1;
for n=3:N
    zfb(n) = exp(1j*(2*pi*fc*T + B*imag(zfb(n-1) - zfb(n-2)))) * zfb(n-1);
end

xfb = zeros(1, N);
xfb(1) = 1;
for n=2:N
    xfb(n) = sin(2*pi*fc*n*T + B*xfb(n-1));
end

Zfb = fft(zfb, Nfft);
ZfbPos = Zfb(1:Nfft/2+1);

Xfb = fft(xfb, Nfft);
XfbPos = Xfb(1:Nfft/2+1);

figure
plot(faxis, 20*log10(abs(ZfbPos)/max(abs(ZfbPos))));
hold on
plot([fc fc], [-60 0], 'r--')
ylim([-60 0]);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Magnitude spectrum of z_{fb}(n)');
set(gca, 'FontSize', 15)
grid on;
legend('z_{fb}', 'f_c')
if saveDissertationFigures==1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, 'figures/FeedbackFMMagSpec', 'epsc')
end


%% Loopback FM (Static Pitch and Timbre)
%%%%%%%%%%%%%%%%%%%

BVec = [1 0.9 0.5 0.0];

for i=1:length(BVec)
    
    B = BVec(i);
    
    zc = zeros(1, N);
    zc(1) = 1;
    for n=2:N
        zc(n) = exp(1j*2*pi*fc*T*(1 + B*real(zc(n-1)))) * zc(n-1);
    end

    f0 = fc * sqrt(1 - B^2);

    Zc = fft(zc, Nfft);
    ZcPos = Zc(1:Nfft/2+1);

    figure
    plot(faxis, 20*log10(abs(ZcPos)/max(abs(ZcPos))));
    hold on
    plot([fc fc], [-60 0], 'r--')
    plot([f0 f0], [-60 0], 'g--')
    xlim([0 1200]);
    ylim([-60 0]);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (dB)');
    title(sprintf('Magnitude spectrum of z_{c}(n), B=%.1f', B));
    set(gca, 'FontSize', 15)
    grid on;
    legend('z_c', 'f_c', 'f_0');
    if saveDissertationFigures==1
        fig = gcf;
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 6 2.5];
        saveas(gcf, ['figures/LoopbackFMMagSpec' num2str(B*10)], 'epsc')
    end
end

% closed form loopback FM
B = 0.99;
fc = 6500;
b0 = (sqrt(1 - B^2) - 1)/B;
f0 = fc*sqrt(1 - B^2);

z0 = (b0 + exp(j*2*pi*f0*nT)) ./ (1 + b0*exp(j*2*pi*f0*nT));

Z0 = fft(z0, Nfft);
Z0Pos = Z0(1:Nfft/2+1);



% loopback FM comparison
zc = zeros(1, N);
zc(1) = 1;
for n=2:N
    zc(n) = exp(j*2*pi*fc*T*(1 + B*real(zc(n-1)))) * zc(n-1);
end

Zc = fft(zc, Nfft);
ZcPos = Zc(1:Nfft/2+1);


figure
subplot(211)
plot(faxis, 20*log10(abs(ZcPos)/max(abs(ZcPos))));
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('z_{c}(n)');
% hold on
% plot([f0 f0], [-60 0], 'r--')
xlim([faxis(1) faxis(end)])
ylim([-60 0])
set(gca, 'FontSize', 12)
grid on;

subplot(212)
plot(faxis, 20*log10(abs(Z0Pos)/max(abs(Z0Pos))));
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('z_{0}(n)');
% hold on
% plot([f0 f0], [-60 0], 'r--')
xlim([faxis(1) faxis(end)])
ylim([-60 0])
set(gca, 'FontSize', 12)
grid on;
sgt = sgtitle(sprintf('Magnitude spectrum comparison for f_c=%d, B=%.2f', fc, B));
sgt.FontSize = 15;
if saveDissertationFigures==1
    fig = gcf;
    saveas(gcf, ['figures/LoopbackFMClosedFormMagSpecComparison'], 'epsc')
end

%% Loopback FM (Time-varying Pitch and Timbre)
%%%%%%%%%%%%%%%%%%%

% pitch glide - parameters
decayT60 = 0.6; % time I want to pitch glide to decay by 60dB
f0_0 = 200;         % starting frequency, Hz
f0_1 = 140;          % ending frequency, Hz

% pitch glide - derived parameters
T60Samp = decayT60*fs;
% exponentially decaying envelope according to the T60
A = 1;
tau = -(decayT60)/log(0.001);
e = A * exp(-nT./tau);
% scaled f0 function
f0Decay = (f0_0 - f0_1) * e + f0_1;

% === loopback FM exponential decay ===
% we need to make sure the w0 <= wc. w0Tilde(1) will always be the largest
% value since this is an exponentially decreasing function. We can use
% w0Tilde(1) and a user-given B to solve for wc.
w0Tilde = 2*pi*f0Decay;
wc = w0Tilde(1)/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0Tilde/wc).^2);

zcTV = zeros(1, N);
zcTV(1) = 1;
for i=2:N
    zcTV(i) = exp(j*wc*T*(1 + BTilde(i) * real(zcTV(i-1)))) * zcTV(i-1);
end

% === Closed form exponential decay ===
% ThetaH = integral [(f0_0 - f0_1) * A * exp(-nT/tau) + f0_1 dnT] <-- that's the equation
b0 = linspace(0.001, 0.999, N);
Theta0 = (2 * pi * (f0_0 - f0_1)) * (-tau * A) * exp(-nT/tau) + (2 * pi * f0_1 * nT); 
z0TV = (b0 + exp(1j.*Theta0)) ./ (1 + b0.*exp(1j.*Theta0));


% === exponential decay pitch glide plots ===
figure

subplot(211)
spectrogram(real(zcTV), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(nT*1000, f0Decay/1000, 'r')
ylim([0 1])
title('z_c(n)');
set(gca, 'FontSize', 12)
subplot(212)
spectrogram(real(z0TV), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(nT*1000, f0Decay/1000, 'r')
ylim([0 1])
title('z_0(n)');
set(gca, 'FontSize', 12)
sgt = sgtitle('Spectrograms of loopback FM with exponentially decaying pitch glides');
sgt.FontSize = 15;
if saveDissertationFigures==1
    fig = gcf;
    saveas(gcf, ['figures/LoopbackFMClosedFormTimeVaryingComparison'], 'epsc')
end
