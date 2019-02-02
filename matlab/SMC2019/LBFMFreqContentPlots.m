% Loopback FM, PM and SAPF frequency content plots

fs = 44100;
dur = 1;

T = 1/fs;
N = fs * dur;
nT = 0:T:(dur-T);

Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

wc = 2*pi*700; 
g = 0.9999;

%% %%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback FM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate magnitude spectrum plot to show how the frequency components are
% integer multiples of w0

h_fm = zeros(1, N);
h_fm(1) = 1;

w0 = 2*pi*300;
B = 0.9;

wc = w0/sqrt(1 - B^2);

for n=2:N
    h_fm(n) = (exp(1j*wc*T*(1 + B * real(h_fm(n-1))))) * h_fm(n-1);
end

H_FM = fft(h_fm, Nfft);
H_FMPos = H_FM(1:Nfft/2+1);
H_FMPosdB = 20*log10(abs(H_FMPos)/max(abs(H_FMPos)));

figure
plot(faxis, H_FMPosdB, 'linewidth', 2)
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Loopback FM magnitude spectrum for \omega _0 = 2\pi300'); 
set(gca, 'FontSize', 15);
ylim([-60 0]);
grid on
saveas(gcf, 'figures/LBFMMagSpec', 'epsc')


%% %%%%%%%%%%%%%%%%%%%%%%%%%
% Stretched Allpass Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BVec = [0.02 0.03 0.04 0.05:0.05:0.95 0.96 0.97 0.98 0.99];
NB = length(BVec);
slopes = zeros(1, NB);
numSidebands = zeros(1, NB);
b0Vec = zeros(1, NB);

for i=1:NB
    
    B = BVec(i);
    w0 = wc*sqrt(1 - B^2);

    Theta = wc*sqrt(1 - B^2)*nT;
    b0 = (sqrt(1 - B^2) - 1)/B;

    H_apf = (b0 + exp(1j.*Theta)) ./ (1 + b0.*exp(1j.*Theta));

    %% get a function for how the peaks decrease as B changes 
    % (what is the slope or rolloff?)



    
    HAPF = fft(H_apf, Nfft);
    HAPFPos = HAPF(1:Nfft/2+1);

    HAPFPosdB = 20*log10(abs(HAPFPos)/max(abs(HAPFPos)));

    % find the peaks dB spectrum plot
    [pks, locs] = findpeaks(HAPFPosdB, 'MinPeakHeight', -60, 'MinPeakDistance', w0/(2*pi));

    % remove first peak at 0 (or some incorrect one that findpeaks() finds
    locs(1) = [];
    pks(1) = [];

    % fit a straight line
    p = polyfit(faxis(locs), pks, 1);

    figure
    plot(faxis, HAPFPosdB, 'linewidth', 2)
    hold on
    for k=1:length(locs)
        plot(faxis(locs(k)), pks(k), 'r*');
    end
    plot(faxis(locs), p(1)*faxis(locs) + p(2), 'g--')
    
    %pause
    close
    
    slopes(i) = p(1);
    numSidebands(i) = length(pks);
    b0Vec(i) = b0;
end

% For SAPF, having this in terms of b0 makes sense
figure
subplot(211)
plot(b0Vec, slopes);
title('b0Vec vs slopes (rolloff)');
xlim([min(b0Vec) max(b0Vec)]);
grid on

subplot(212)
plot(b0Vec, numSidebands);
title('b0Vec vs number sidebands');
xlim([min(b0Vec) max(b0Vec)]);
grid on

% fit that to an exponential function
b00 = -0.8676; y00 = 40;
b01 = -0.01; y01 = 2;
A0 = ((y01^(b00/b01))/y00)^(1/(b00/b01-1));
tau = -b01/(log(y01/A0));
y = A0*exp(-b0Vec/tau);

hold on
plot(b0Vec, y);
legend('data', 'fit')

% for LBFM and LBPM, having this in terms of B makes sens
figure
subplot(211)
plot(BVec, slopes);
title('BVec vs slopes (rolloff)');
xlim([min(BVec) max(BVec)]);
grid on

% fit that to a logarithmic growth function
x0 = 0.02; y0 = -0.05783;
x1 = 0.99; y1 = -0.01253;
a = (y1*log(x0) - y0*log(x1))/(log(x0) - log(x1));
b = (y0 - a)/log(x0);
y = a + b*log(BVec);
hold on
plot(BVec, y);
legend('data', 'fit')

subplot(212)
plot(BVec, numSidebands);
title('BVec vs number sidebands');
xlim([min(BVec) max(BVec)]);
grid on

% fit that to a double exponential growth function (not the best fit)
x0 = 0.02; y0 = 2;
x1 = 0.99; y1 = 40;
k = 6.6;
y = y0 + (exp(k.^BVec)/max(exp(k.^BVec)))*(y1-y0);

hold on
plot(BVec, y);


%% %%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
subplot(311)
spectrogram(H_apf, hann(1024), 512, 1024, fs, 'yaxis')
hold on
plot(linspace(0, 1000, N), w0/(2*pi*1000) * ones(1, N), 'r')
ylim([0 2])
title(sprintf('Allpass Filter with B=%f and f0=%f', B, w0/(2*pi)));
legend('w_0')

subplot(312)
spectrogram(h_fm, hann(1024), 512, 1024, fs, 'yaxis')
hold on
plot(linspace(0, 1000, N), w0/(2*pi*1000) * ones(1, N), 'r')
ylim([0 2])
title(sprintf('Feedback FM with B=%f and f0=%f', B, w0/(2*pi)));
legend('w_0')

subplot(313)
spectrogram(h_pm, hann(1024), 512, 1024, fs, 'yaxis')
hold on
plot(linspace(0, 1000, N), w0/(2*pi*1000) * ones(1, N), 'r') 
ylim([0 2])
title(sprintf('Feedback PM B=%f and f0=%f', B, w0/(2*pi)));
legend('w_0')

