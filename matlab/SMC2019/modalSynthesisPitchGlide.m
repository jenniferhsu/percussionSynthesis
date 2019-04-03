% modalSynthesisTests.m
% Modal Synthesis as implemented as bandpass filters and with pitch glides
% If we are doing modal synthesis with bandpass filters, the pitch glides
% are straightforward. Use whatever frequency glide you want (set f0 to the
% frequency glide) and directly set f_k in the bandpass filters to the f0
% function.  When we do pitch glides for sinusoids, it's more annoying
% because we have to take the integral of the phase to get the correct
% sounding frequencies.

% parameters
fs = 44100;
dur = 1;
f0 = 440;   % frequency glide (first frequency)
f1 = 220;   % frequency glide (second frequency)
r = [1];    % decays
K = 1;      % number of modes

% derived parameters
N = fs*dur;
T = 1/fs;
nT = 0:T:(dur-T);

% linear pitch glide
m = f1 - f0;
b = (f1 - m*dur);
fL = m * nT + b;
ffL = m/2*nT.^2 + b*nT;  

% exponential pitch glide
T60 = 0.8;
A = 1;
tau = -T60/log(0.001);
e = A*exp(-nT/tau);
fE = abs(f1 - f0)*e + min(f0, f1);
ffE = -abs(f1 - f0) * tau * A * exp(-nT/tau) + min(f0, f1)*nT;

% fft parameters
Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

x = zeros(N, 1);
x(1) = 1;

% bandpass filter with linear pitch glide
zL = [0 0];
yL = zeros(1, N);
for k=1:K
    B = 1;
    for n=1:N
        %A = [1 -2*r(k)*cos(2*pi*f(k)*T) r(k)^2];
        A = [1 -2*r(k)*cos(2*pi*fL(k,n)*T) r(k)^2];
        [yL(n), zL] = filter(B, A, x(n), zL);
    end
end

zL2 = [0 0];
yL2 = zeros(1, N);
for k=1:K
    B = 1;
    for n=1:N
        %A = [1 -2*r(k)*cos(2*pi*f(k)*T) r(k)^2];
        A = [1 -2*r(k)*cos(2*pi*ffL(k,n)*T) r(k)^2];
        [yL2(n), zL2] = filter(B, A, x(n), zL2);
    end
end

figure
subplot(211)
spectrogram(real(yL), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1, N)*1000, fL/1000, 'r');
ylim([0 1])
title('bandpass filter with linear pitch glide (correct)')
subplot(212)
spectrogram(real(yL2), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1, N)*1000, fE/1000, 'r');
ylim([0 1])
title('bandpass filter with linear pitch glide (using the integral)')


% bandpass filter with exponential pitch glide
zE = [0 0];
yE = zeros(1, N);
for k=1:K
    B = 1;
    for n=1:N
        %A = [1 -2*r(k)*cos(2*pi*f(k)*T) r(k)^2];
        A = [1 -2*r(k)*cos(2*pi*fE(k,n)*T) r(k)^2];
        [yE(n), zE] = filter(B, A, x(n), zE);
    end
end

zE2 = [0 0];
yE2 = zeros(1, N);
for k=1:K
    B = 1;
    for n=1:N
        %A = [1 -2*r(k)*cos(2*pi*f(k)*T) r(k)^2];
        A = [1 -2*r(k)*cos(2*pi*ffE(k,n)*T) r(k)^2];
        [yE2(n), zE2] = filter(B, A, x(n), zE2);
    end
end


figure
subplot(211)
spectrogram(real(yE), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1, N)*1000, fE/1000, 'r');
ylim([0 1])
title('bandpass filter with linear pitch glide (correct)')
subplot(212)
spectrogram(real(yE2), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1, N)*1000, fE/1000, 'r');
ylim([0 1])
title('bandpass filter with exponential pitch glide')




% check with a sinusoid
%y2 = sin(2*pi*f.*nT);
%y2 = exp(1j*2*pi*f.*nT);

% linear pitch glide integral
y2 = sin(2*pi*ff);
%y2 = exp(1j*2*pi*ff);

% exponential pitch glide integral
y21 = sin(2*pi*f);
y22 = sin(2*pi*ff);
%y2 = exp(1j*2*pi*ff);

figure
subplot(211)
spectrogram(real(y21), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1, N)*1000, f/1000, 'r');
ylim([0 1])
title('sine oscillator with pitch glide (incorrect)')
subplot(212)
spectrogram(real(y22), hann(256), 128, 1024, fs, 'yaxis');
hold on
plot(linspace(0, 1, N)*1000, f/1000, 'r');
ylim([0 1])
title('sine oscillator with pitch glide (correct)')