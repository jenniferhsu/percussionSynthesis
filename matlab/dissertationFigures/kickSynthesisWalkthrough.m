% Kick example
saveDissertationFigures = 1;
f0_0 = 100;
f0_1 = 40;
r = 0.6;

fs = 44100;
dur = 1;

N = fs*dur;
T = 1/fs;
nT = [0:1:N-1] * T;

% solve for exponentially decaying pitch glide envelope
A = 1;
tau = -r/(log(0.001));
d = exp(-nT/tau);

%% get the pitch contour
f0Tilde = (f0_0 - f0_1)*d + f0_1;

figure
plot(nT, f0Tilde, 'linewidth', 2);
xlabel('Time (seconds)');
ylabel('Amplitude (linear)');
title('Kick Drum Example: $\tilde{f}_0(n)$','Interpreter','latex')
grid on
set(gca, 'FontSize', 15)
if saveDissertationFigures==1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, 'figures/KickDrumf0Tilde', 'epsc')
end

%% with loopback FM oscillators (z_c(n))
w0Tilde = 2 * pi * f0Tilde;
wc = w0Tilde(1);
BTilde = sqrt(1 - (w0Tilde / wc).^2);

zc = zeros(1, N);
zc(1) = 1;
for i=2:N
    zc(i) = exp(j*wc*T*(1 + BTilde(i) * real(zc(i-1)))) * zc(i-1);
end

% plot BTilde
figure
plot(nT, BTilde, 'linewidth', 2);
xlabel('Time (seconds)');
ylabel('Amplitude (linear)');
title('Kick Drum Example: $\tilde{B}(n)$','Interpreter','latex')
grid on
set(gca, 'FontSize', 15)
if saveDissertationFigures==1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, 'figures/KickDrumBTilde', 'epsc')
end

%% with closed form loopback FM (z_0(n))
Theta0 = -2 * pi * (f0_0 - f0_1) * tau * A* exp(-nT/tau) + 2 * pi * f0_1 * nT;
%b0 = (sqrt(1 - BTilde.^2) - 1) ./ BTilde;
b0 = 0.1;
z0 = (b0 + exp(1j.*Theta0)) ./ (1 + b0.*exp(1j.*Theta0));

% plot pitch glide with z_0(n) and z_c(n)
figure
subplot(211)
spectrogram(real(zc), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(nT*1000, f0Tilde/1000, 'r')
ylim([0 1])
title('z_c(n)');
set(gca, 'FontSize', 12)
colorbar('off')
subplot(212)
spectrogram(real(z0), hann(1024), 512, 2048, fs, 'yaxis')
hold on
plot(nT*1000, f0Tilde/1000, 'r')
ylim([0 1])
title('z_0(n) with b_0 = 0.1');
set(gca, 'FontSize', 12)
colorbar('off')
sgt = sgtitle('Kick Drum Example: Spectrogram of pitch glide with z_c(n) and z_0(n)');
sgt.FontSize = 15;
if saveDissertationFigures==1
    fig = gcf;
    saveas(gcf, ['figures/KickDrumPitchGlide'], 'epsc')
end

%% set a decay time
T60 = 0.7;
A = 1;

lambda = -T60/log(0.001);

