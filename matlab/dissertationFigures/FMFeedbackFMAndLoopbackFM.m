fs = 44100;
dur = 1;
saveDissertationFigures = 1;

% carrier and modulator frequencies
fc = 200;
fm = 150;
I = 1;

% derived parameters
N = fs*dur;
T = 1/fs;

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
title('Magnitude spectrum of z_{fmc}(n)');
set(gca, 'FontSize', 15)
grid on;
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
B = 0.9;
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
title('Magnitude spectrum of z_{fmc}(n)');
set(gca, 'FontSize', 15)
grid on;
if saveDissertationFigures==1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, 'figures/FeedbackFMMagSpec', 'epsc')
end