addpath(genpath('../proofOfConcept'));

savePlots = 1;
plotOutDir = 'figures/';

% input parameters
fs = 44100;
dur = 0.9;
decayT60 = 0.6;     % time i want the pitch glide envelope to be at -60dB
f0 = 1000;

% derived parameters
N = dur*fs;
T = 1/fs;
n = 0:1:N-1;
nT = 0:T:(dur-T);
w0 = 2*pi*f0;

% decay envelope
A_e = 1;
tau_e = -(N-1)/log(0.001);
env = A_e * exp(-n/tau_e);


%% CLOSED FORM LOOPBACK FM LINEARLY INCREASING b0

b00 = 0;
b01 = 0.8;
k = (b01 - b00)/N;
l = b01 - k*N;
b0 = k * n + l;
y0Lin = (b0 + exp(1j*w0*n*T)) ./ (1 + b0.*exp(1j*w0*n*T));


%% LOOPBACK FM LINEARLY INCREASING b0 so B is different

BTilde = - 2*b0 ./(b0.^2 + 1);
ycLin = zeros(1, N);
ycLin(1) = 1;
for i=2:N
    ycLin(i) = exp(1j*wc*T*(1 + BTilde(i) * real(ycLin(i-1)))) * ycLin(i-1);
end


%% linearly increasing b0(n) plots

figure
subplot(211)
spectrogram(real(y0Lin), hann(1024), 512, 2048, fs, 'yaxis')
%ylim([0 1])
title('Timbre variation with b_0(n) increasing linearly with z_0(n)');
colorbar('off')
subplot(212)
spectrogram(real(ycLin), hann(1024), 512, 2048, fs, 'yaxis')
%ylim([0 1])
title('Timbre variation with b_0(n) increasing linearly with z_0(n)');
colorbar('off')

if savePlots
    H = figure
    subplot(211)
    spectrogram(real(y0Lin), hann(1024), 512, 2048, fs, 'yaxis')
    title('Timbre variation with linearly increasing $b_0(n)$ using z_0(n)');
    set(gca, 'fontsize', 15)
    subplot(212)
    spectrogram(real(ycLin), hann(1024), 512, 2048, fs, 'yaxis')
    title('Timbre variation with linearly increasing $b_0(n)$ using z_c(n)');
    set(gca, 'fontsize', 15)
    saveas(gcf, [plotOutDir 'timbreVariationLinIncrease'], 'epsc')
%     fig = H
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 6 2.6];
%     print([plotOutDir 'timbreVariationLinIncrease'], '-depsc', '-r0')

    figure
    subplot(211)
    plot(nT, b0, 'linewidth', 2)
    grid on
    xlabel('Time (seconds)')
    ylabel('b_0(n)')
    subplot(212)
    plot(nT, BTilde, 'linewidth', 2)
    grid on
    xlabel('Time (seconds)')
    ylabel('B(n)')
    sgtitle('Linearly increasing b_0(n) and corresponding B(n) function');
    saveas(gcf, [plotOutDir 'timbreVariationLinb0B'], 'epsc')
end

%% CLOSED FORM LOOPBACK FM EXPONENTIALLY DECREASING b0

b00 = 0;
b01 = 0.8;
b0 = (b01 - b00) * env + b00;

y0Exp = (b0 + exp(1j*w0*n*T)) ./ (1 + b0.*exp(1j*w0*n*T));


%% LOOPBACK FM EXPONENTIALLY DECREASING b0 so B is different

BTilde = - 2*b0 ./(b0.^2 + 1);
ycExp = zeros(1, N);
ycExp(1) = 1;
for i=2:N
    ycExp(i) = exp(1j*wc*T*(1 + BTilde(i) * real(ycExp(i-1)))) * ycExp(i-1);
end

%% exponentially decreasing b0(n) plots

figure
subplot(211)
spectrogram(real(y0Exp), hann(1024), 512, 2048, fs, 'yaxis')
%ylim([0 1])
title('Timbre variation with exponentially decreasing b_0(n) with z_0(n)');
colorbar('off')
subplot(212)
spectrogram(real(ycExp), hann(1024), 512, 2048, fs, 'yaxis')
%ylim([0 1])
title('Timbre variation with exponentially decreasing b_0(n) with z_c(n)');
colorbar('off')

if savePlots
    H = figure
    subplot(211)
    spectrogram(real(y0Exp), hann(1024), 512, 2048, fs, 'yaxis')
    title('Timbre variation with exponentially decreasing $b_0(n)$ using z_0(n)');
    set(gca, 'fontsize', 15)
    subplot(212)
    spectrogram(real(ycExp), hann(1024), 512, 2048, fs, 'yaxis')
    title('Timbre variation with exponentially decreasing $b_0(n)$ using z_c(n)');
    set(gca, 'fontsize', 15)
    saveas(gcf, [plotOutDir 'timbreVariationExpDecrease'], 'epsc')
%     fig = H
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 6 2.6];
%     print([plotOutDir 'timbreVariationLinIncrease'], '-depsc', '-r0')

    figure
    subplot(211)
    plot(nT, b0, 'linewidth', 2)
    grid on
    xlabel('Time (seconds)')
    ylabel('b_0(n)')
    subplot(212)
    plot(nT, BTilde, 'linewidth', 2)
    grid on
    xlabel('Time (seconds)')
    ylabel('B(n)')
    sgtitle('Exponentially decreasing b_0(n) and corresponding B(n) function');
    saveas(gcf, [plotOutDir 'timbreVariationExpb0B'], 'epsc')
end
