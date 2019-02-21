fs = 44100;
dur = 1;
wc = 2*pi*700;
B = 0.9;

N = fs*dur;
T = 1/fs;
nT = linspace(0, dur, N);

g = 0.9999;
BVec = g.^(0:N-1)';
w0Vec = wc*sqrt(1 - BVec.^2);
b0Vec = (sqrt(1 - BVec.^2) - 1)./BVec;

% Loopback FM
h_fm = zeros(1, N);
h_fm(1) = 1;
for n=2:N
    h_fm(n) = exp(1j*wc*T*(1 + BVec(n).*real(h_fm(n-1)))) * h_fm(n-1);
end

% Loopback PM
h_pm = zeros(1, N);
h_pm(1) = 1;
s = 0;
for n=2:N
    s = s + (BVec(n-1) * h_pm(n-1));
    h_pm(n) = exp(1j*wc*T*(n + real(s)));
end

% stretched APF
C = 0;
u = sqrt(1 - BVec.^2);
Theta = (wc*T/log(g)) .* (u - atanh(u) + C);

H = (b0Vec + exp(1j*Theta)) ./ (1 + b0Vec .* exp(1j*Theta));

% stretched APF, angle  version
%angleH = w0Vec.*nT' - 2*atan((b0Vec .* sin(w0Vec.*nT')) ./ (1 + b0Vec .* cos(w0Vec.*nT')));
% above is wrong
angleH = angle(H);
H1 = exp(1j*angleH);

% Incorrect Loopback PM
h_pm1 = zeros(1, N);
h_pm1(1) = 1;
s = 0;
for n=2:N
    s = s + h_pm1(n-1);
    h_pm1(n) = exp(1j*wc*T*(n + BVec(n-1)*real(s)));
end

figure

subplot(311)
spectrogram(h_fm, hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1000, N), w0Vec/(2*pi*1000), 'r')
title('time-varying loopback FM');
ylim([0 10])

% subplot(312);
% spectrogram(h_pm, hann(256), 128, 1024, fs, 'yaxis');
% hold on
% plot(linspace(0, 1000, N), w0Vec/(2*pi*1000), 'r')
% title('time-varying loopback PM');
% ylim([0 10])

subplot(312)
spectrogram(H, hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1000, N), w0Vec/(2*pi*1000), 'r')
title('time-varying stretched APF');
ylim([0 10])

subplot(313)
spectrogram(H1, hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1000, N), w0Vec/(2*pi*1000), 'r')
title('time-varying stretched APF');
ylim([0 10])


figure
subplot(211)
spectrogram(h_pm1, hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1000, N), w0Vec/(2*pi*1000), 'r')
title('incorrect time-varying loopback PM');
ylim([0 10])

subplot(212)
spectrogram(H1, hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1000, N), w0Vec/(2*pi*1000), 'r')
title('incorrect time-varying stretched APF, angle version');
ylim([0 10])
