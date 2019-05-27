% pitchGlideB.m
% linear and exponential B(n) plots


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


%% envelope
fadeTime = 0.1; % seconds
fadeSamps = fadeTime*fs;
env = [ones(1, N-fadeSamps) linspace(1, 0, fadeSamps)];

%% Linear increase for zc(n)

B0 = 0;
B1 = 0.9999;
k = (B1 - B0)/N;
l = 0;

% zc(n)
wc = 2*pi*840;
BTilde = k*n + l;

yLin = zeros(1, N);
yLin(1) = 1;
for i=2:N
    yLin(i) = exp(j*wc*T*(1 + BTilde(i) * real(yLin(i-1)))) * yLin(i-1);
end

%% find the pitch glide
w0Lin = wc * sqrt(1 - BTilde.^2);
f0Lin = w0Lin/(2*pi);

%% Linear increase for z0(n)

b0 = -0.3;
Theta0 = wc*T/(2*k)*(BTilde .* sqrt(1 - BTilde.^2) + asin(BTilde));
YLin = (b0 + exp(j*Theta0)) ./ (1 + b0*exp(j*Theta0));

%% linear B(n) pitch glide plots


if savePlots
    
    figure
    subplot(211)
    plot(n*T*1000, real(Theta0), 'linewidth', 2)
    grid on
    xlabel('Time (ms)');
    ylabel('\Theta_0(n) value');
    title('\Theta_0(n)');
    set(gca, 'fontsize', 15)
    
    subplot(212)
    plot(n*T, BTilde, 'linewidth', 2)
    grid on
    xlabel('Time (ms)');
    ylabel('B(n) value');
    title('B(n)');
    set(gca, 'fontsize', 15)
    saveas(gcf, [plotOutDir 'pitchGlideLinIncrease_forB_BAndTheta0'], 'epsc')
    
    H = figure
    subplot(211)
    spectrogram(real(YLin), hann(1024), 512, 2048, fs, 'yaxis')
    colorbar('off')
    hold on
    plot(n*T*1000, f0Lin/1000, 'r', 'linewidth', 2)
    ylim([0 10])
    title('Pitch glide increasing B(n) linearly with z_0(n)');
    set(gca, 'fontsize', 15)
    subplot(212)
    spectrogram(real(yLin), hann(1024), 512, 2048, fs, 'yaxis')
    colorbar('off')
    hold on
    plot(n*T*1000, f0Lin/1000, 'r', 'linewidth', 2)
    ylim([0 10])
    title('Pitch glide increasing B(n) linearly with z_c(n)');
    set(gca, 'fontsize', 15)
    saveas(gcf, [plotOutDir 'pitchGlideLinBIncrease'], 'epsc')
%     fig = H
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 6 2.6];
%     print([plotOutDir 'pitchGlideLinIncrease'], '-depsc', '-r0')
    %close
end


%% save audio
audiowrite([audioDir 'pitchGlideLinIncreaseB_z0.wav'], scaleForSavingAudio(real(YLin) .* env), fs);
audiowrite([audioDir 'pitchGlideLinIncreaseB_zc.wav'], scaleForSavingAudio(real(yLin) .* env), fs);

%% Exponential increase for zc(n)

g = 0.9999;

% zc(n)
wc = 2*pi*840;
BTilde = g.^n;

yExp = zeros(1, N);
yExp(1) = 1;
for i=2:N
    yExp(i) = exp(j*wc*T*(1 + BTilde(i) * real(yExp(i-1)))) * yExp(i-1);
end

%% find the pitch glide
w0Exp = wc * sqrt(1 - BTilde.^2);
f0Exp = w0Exp/(2*pi);

%% Exponential increase for z0(n)

b0 = -0.3;
Theta0 = (wc*T/log(g))*(sqrt(1-g.^(2*n)) - atanh(sqrt(1-g.^(2*n))));
YExp = (b0 + exp(j*Theta0)) ./ (1 + b0*exp(j*Theta0));


%% Exponential pitch glide plots

if savePlots
    
      figure
    subplot(211)
    plot(n*T*1000, real(Theta0), 'linewidth', 2)
    grid on
    xlabel('Time (ms)');
    ylabel('\Theta_0(n) value');
    title('\Theta_0(n)');
    set(gca, 'fontsize', 15)
    
    subplot(212)
    plot(n*T, BTilde, 'linewidth', 2)
    grid on
    xlabel('Time (ms)');
    ylabel('B(n) value');
    title('B(n)');
    set(gca, 'fontsize', 15)
    saveas(gcf, [plotOutDir 'pitchGlideExpDecay_forB_BAndTheta0'], 'epsc')
    
    H = figure
    subplot(211)
    spectrogram(real(YExp), hann(1024), 512, 2048, fs, 'yaxis')
    hold on
    plot(n*T*1000, f0Exp/1000, 'r', 'linewidth', 2)
    ylim([0 10])
    colorbar('off')
    title('Pitch glide for exponentially decreasing B(n) with z_0(n)');
    set(gca, 'fontsize', 15)
    subplot(212)
    spectrogram(real(yExp), hann(1024), 512, 2048, fs, 'yaxis')
    hold on
    plot(n*T*1000, f0Exp/1000, 'r', 'linewidth', 2)
    ylim([0 10])
    colorbar('off')
    title('Pitch glide for exponentially decreasing B(n) with z_c(n)');
    set(gca, 'fontsize', 15)
    saveas(gcf, [plotOutDir 'pitchGlideExpDecayB'], 'epsc')
%     fig = H
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 6 2.6];
%     print([plotOutDir 'pitchGlideExpDecay'], '-depsc', '-r0')
    %close
end



%% save audio
audiowrite([audioDir 'pitchGlideExpDecreaseB_z0.wav'], scaleForSavingAudio(real(YExp) .* env), fs);
audiowrite([audioDir 'pitchGlideExpDecreaseB_zc.wav'], scaleForSavingAudio(real(yExp) .* env), fs);