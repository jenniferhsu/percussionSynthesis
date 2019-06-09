% pitchGlidef0.m
% linear and exponential f0(n) plots

addpath(genpath('../proofOfConcept'));

savePlots = 1;
plotOutDir = 'figures/';
audioDir = 'audioExamples/';

% input parameters
fs = 44100;
dur = 0.9;
decayT60 = 0.6;     % time i want the pitch glide envelope to be at -60dB
b0 = -0.3;        % timbre control
B = -2*b0 /(b0^2 + 1); 

% derived parameters
N = dur*fs;
T = 1/fs;
n = 0:1:N-1;
nT = 0:T:(dur-T);
T60Samp = decayT60*fs;

% decay envelope
A_e = 1;
tau_e = -(N-1)/log(0.001);
env = A_e * exp(-n/tau_e);

%w0 = 2*pi*f0_0;

col = colormap(pink)

%% envelope
fadeTime = 0.1; % seconds
fadeSamps = fadeTime*fs;
env = [ones(1, N-fadeSamps) linspace(1, 0, fadeSamps)];


%% Linear increase

f0_0 = 140;         % starting frequency, Hz
f0_1 = 840;          % ending frequency, Hz

m = (f0_1 - f0_0)/dur;
b = f0_0;
f0Lin = m*nT + b;
w0Tilde = 2*pi*f0Lin;

%% LOOPBACK FM LINEAR INCREASE

wc = w0Tilde(end)/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0Tilde/wc).^2);

yLin = zeros(1, N);
yLin(1) = 1;
for i=2:N
    yLin(i) = exp(j*wc*T*(1 + BTilde(i) * real(yLin(i-1)))) * yLin(i-1);
end

%% CLOSED FORM LOOPBACK FM LINEAR INCREASE - pitch glide only

b0 = -0.5;
% Theta0 = integral[2*pi * (m*nT + b) dnT]
Theta0 = 2*pi * (m/2 * nT.^2 + b*nT);
YLinPitchGlide = (b0 + exp(1j.*Theta0)) ./ (1 + b0.*exp(1j.*Theta0));

%% CLOSED FORM LOOPBACK FM LINEAR INCREASE - harmonics change only

b0 = linspace(0.01, -0.999, N);
% Theta0 = integral[2*pi * (m*nT + b) dnT]
Theta0 = 2*pi * (m/2 * nT.^2 + b*nT);
YLinHarmonicsChange = (b0 + exp(1j.*Theta0)) ./ (1 + b0.*exp(1j.*Theta0));


%% linearly increasing pitch glide plots


if savePlots
    
%     figure
%     subplot(211)
%     plot(n*T*1000, real(Theta0), 'linewidth', 2)
%     grid on
%     xlabel('Time (ms)');
%     ylabel('\Theta_0(n) value');
%     title('\Theta_0(n)');
%     set(gca, 'fontsize', 15)
%     
%     subplot(212)
%     plot(n*T, BTilde, 'linewidth', 2)
%     grid on
%     xlabel('Time (ms)');
%     ylabel('B(n) value');
%     title('B(n)');
%     set(gca, 'fontsize', 15)
%     saveas(gcf, [plotOutDir 'pitchGlideLinIncreaseBAndTheta0'], 'epsc')
    
    figure
    spectrogram(real(YLinPitchGlide), hann(1024), 512, 2048, fs, 'yaxis')
    colormap('pink');
    colorbar('off')
    hold on
    plot(n*T*1000, f0Lin/1000, 'color', [0.5 0.5 0.5], 'linewidth', 2)
    ylim([0 8])
    %title('Pitch glide increasing linearly with z_0(n)');
    set(gca, 'fontsize', 13)
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    print([plotOutDir 'pitchGlideYLinPitchGlide'], '-depsc', '-r0')
    
    figure
    spectrogram(real(YLinHarmonicsChange), hann(1024), 512, 2048, fs, 'yaxis')
    colormap('pink');
    colorbar('off')
    hold on
    plot(n*T*1000, f0Lin/1000, 'color', [0.5 0.5 0.5], 'linewidth', 2)
    ylim([0 8])
    % title('Pitch glide increasing linearly with z_c(n)');
    set(gca, 'fontsize', 15)
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    print([plotOutDir 'pitchGlideYLinHarmonicsChange'], '-depsc', '-r0')
    
    figure
    spectrogram(real(yLin), hann(1024), 512, 2048, fs, 'yaxis')
    colormap('pink');
    colorbar('off')
    hold on
    plot(n*T*1000, f0Lin/1000, 'color', [0.5 0.5 0.5], 'linewidth', 2)
    ylim([0 8])
    % title('Pitch glide increasing linearly with z_c(n)');
    set(gca, 'fontsize', 15)
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    print([plotOutDir 'pitchGlideyLin'], '-depsc', '-r0')

end


%% save audio
audiowrite([audioDir 'pitchGlideLinIncrease_z0HarmonicsChange.wav'], scaleForSavingAudio(real(YLinHarmonicsChange) .* env), fs);
audiowrite([audioDir 'pitchGlideLinIncrease_z0PitchGlide.wav'], scaleForSavingAudio(real(YLinPitchGlide) .* env), fs);
audiowrite([audioDir 'pitchGlideLinIncrease_zc.wav'], scaleForSavingAudio(real(yLin) .* env), fs);


%% EXPONENTIAL DECAY

f0_0 = 840;         % starting frequency, Hz
f0_1 = 140;          % ending frequency, Hz

% exponentially decaying envelope according to the T60
A = 1;
tau = -(decayT60)/log(0.001);
e = A * exp(-nT./tau);

% scaled f0 function
f0Exp = (f0_0 - f0_1) * e + f0_1;

figure
plot(nT, f0Exp);
hold on
plot([decayT60 decayT60], [f0_0 f0_1], 'k--')
plot([0 decayT60], [f0_1 f0_1], 'k--')


%% LOOPBACK FM EXPONENTIAL DECAY

% loopback FM variables
% we need to make sure the w0 <= wc. w0Tilde(1) will always be the largest
% value since this is an exponentially decreasing function. We can use
% w0Tilde(1) and a user-given B to solve for wc.
w0Tilde = 2*pi*f0Exp;
wc = w0Tilde(1)/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0Tilde/wc).^2);

yExp = zeros(1, N);
yExp(1) = 1;
for i=2:N
    yExp(i) = exp(j*wc*T*(1 + BTilde(i) * real(yExp(i-1)))) * yExp(i-1);
end

%soundsc(real(yDec) .* env, fs)

%% CLOSED FORM LOOPBACK FM EXPONENTIAL DECAY

b0 = -0.5;
% Theta0 = integral [(f0_0 - f0_1) * A * exp(-nT/tau) + f0_1 dnT] <-- that's the equation
Theta0 = (2 * pi * (f0_0 - f0_1)) * (-tau * A) * exp(-nT/tau) + (2 * pi * f0_1 * nT); 

%b0 = (sqrt(1 - B^2) - 1)/B;
YExpPitchGlide = (b0 + exp(1j.*Theta0)) ./ (1 + b0.*exp(1j.*Theta0));

%% CLOSED FORM LOOPBACK FM HARMONICS ONLY

b0 = linspace(0.01, -0.999, N);
% Theta0 = integral [(f0_0 - f0_1) * A * exp(-nT/tau) + f0_1 dnT] <-- that's the equation
Theta0 = (2 * pi * (f0_0 - f0_1)) * (-tau * A) * exp(-nT/tau) + (2 * pi * f0_1 * nT); 

%b0 = (sqrt(1 - B^2) - 1)/B;
YExpHarmonicsChange = (b0 + exp(1j.*Theta0)) ./ (1 + b0.*exp(1j.*Theta0));

%% decaying pitch glide plots
% 
% 
if savePlots
%     
%     figure
%     subplot(211)
%     plot(n*T*1000, real(Theta0), 'linewidth', 2)
%     grid on
%     xlabel('Time (ms)');
%     ylabel('\Theta_0(n) value');
%     title('\Theta_0(n)');
%     set(gca, 'fontsize', 15)
%     
%     subplot(212)
%     plot(n*T, BTilde, 'linewidth', 2)
%     grid on
%     xlabel('Time (ms)');
%     ylabel('B(n) value');
%     title('B(n)');
%     set(gca, 'fontsize', 15)
%     saveas(gcf, [plotOutDir 'pitchGlideExpDecayBAndTheta0'], 'epsc')
%     
%     H = figure
%     subplot(211)
%     spectrogram(real(YExp), hann(1024), 512, 2048, fs, 'yaxis')
%     colorbar('off')
%     hold on
%     plot(n*T*1000, f0Exp/1000, 'r', 'linewidth', 2)
%     ylim([0 2])
%     colorbar('off')
%     title('Pitch glide decreasing exponentially with z_0(n)');
%     set(gca, 'fontsize', 15)
%     subplot(212)
%     spectrogram(real(yExp), hann(1024), 512, 2048, fs, 'yaxis')
%     colorbar('off')
%     hold on
%     plot(n*T*1000, f0Exp/1000, 'r', 'linewidth', 2)
%     ylim([0 2])
%     colorbar('off')
%     title('Pitch glide decreasing exponentially with z_c(n)');
%     set(gca, 'fontsize', 15)
%     saveas(gcf, [plotOutDir 'pitchGlideExpDecay'], 'epsc')

    figure
    spectrogram(real(YExpPitchGlide), hann(1024), 512, 2048, fs, 'yaxis')
    colormap('pink');
    colorbar('off')
    hold on
    plot(n*T*1000, f0Exp/1000, 'color', [0.5 0.5 0.5], 'linewidth', 2)
    ylim([0 8])
    %title('Pitch glide increasing linearly with z_0(n)');
    set(gca, 'fontsize', 13)
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    print([plotOutDir 'pitchGlideYExpPitchGlide'], '-depsc', '-r0')
    
    figure
    spectrogram(real(YExpHarmonicsChange), hann(1024), 512, 2048, fs, 'yaxis')
    colormap('pink');
    colorbar('off')
    hold on
    plot(n*T*1000, f0Exp/1000, 'color', [0.5 0.5 0.5], 'linewidth', 2)
    ylim([0 8])
    % title('Pitch glide increasing linearly with z_c(n)');
    set(gca, 'fontsize', 15)
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    print([plotOutDir 'pitchGlideYExpHarmonicsChange'], '-depsc', '-r0')
    
    figure
    spectrogram(real(yExp), hann(1024), 512, 2048, fs, 'yaxis')
    colormap('pink');
    colorbar('off')
    hold on
    plot(n*T*1000, f0Exp/1000, 'color', [0.5 0.5 0.5], 'linewidth', 2)
    ylim([0 8])
    % title('Pitch glide increasing linearly with z_c(n)');
    set(gca, 'fontsize', 15)
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    print([plotOutDir 'pitchGlideyExp'], '-depsc', '-r0')

end

%% save audio
audiowrite([audioDir 'pitchGlideExpDecay_z0HarmonicsChange.wav'], scaleForSavingAudio(real(YExpHarmonicsChange) .* env), fs);
audiowrite([audioDir 'pitchGlideExpDecay_z0PitchGlide.wav'], scaleForSavingAudio(real(YExpPitchGlide) .* env), fs);
audiowrite([audioDir 'pitchGlideExpDecay_zc.wav'], scaleForSavingAudio(real(yExp) .* env), fs);