

% parameters
fs = 44100;
dur = 1;
fc = 220;       % carrier frequency, Hz
B = 0.5;        % modulation index

% derived parameters
N = fs*dur;
T = 1/fs;
nT = 0:T:(dur-T);

% fft parameters
Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

% feedback FM (discrete formula from Andrew Horner's 1998 paper called
% "Nested Modulator and Feedback FM Matching of Instrument Tones)
x = zeros(1, N);
for n=0:N-1
    ind = n+1;
    if ind-1 < 1
        x(ind) = sin((2*pi*fc*n/fs) + 0);
    else
        x(ind) = sin((2*pi*fc*n/fs) + B*x(ind-1));
    end
end

% feedback FM calculated similarly to loopback FM
x2 = zeros(1, N);
x2(2) = 1;
for n=2:N-1
    ind = n + 1;
    x2(ind) = exp(1j*2*pi*fc*T) * exp(1j*B*imag(x2(ind-1) - x2(ind-2))) * x2(ind-1);
end

% loopback FM as a stretched allpass filter
w0 = 2*pi*fc;
wc = w0/sqrt(1 - B^2);
b0 = (sqrt(1 - B^2) - 1)/B;
h = (b0 + exp(1j*w0*nT))./(1 + b0*exp(1j*w0*nT));

% loopback FM as a rotation
h = zeros(1, N);
h(1) = 1;
for n=1:N-1
    ind = n+1;
    h(ind) = exp(1j*wc*T*(1 + B*real(h(ind-1)))) * h(ind-1);
end

% fft
%X = fft(x, Nfft);
X = fft(x2, Nfft);
XPos = X(1:Nfft/2+1);
XPosdB = 20*log10(abs(XPos)/max(abs(XPos)));

H = fft(real(h), Nfft);
HPos = H(1:Nfft/2+1);
HPosdB = 20*log10(abs(HPos)/max(abs(HPos)));


figure
subplot(211)
plot(faxis, XPosdB);
ylim([-60 0])
xlim([0 1000])
%spectrogram(real(x), hann(256), 128, 1024, fs, 'yaxis');
%title('feedback FM spectrogram')
%ylim([0 1])
subplot(212)
plot(faxis, HPosdB);
ylim([-60 0])
xlim([0 1000])
%spectrogram(real(h), hann(256), 128, 1024, fs, 'yaxis');
%title('loopback FM spectrogram')
%ylim([0 1])
