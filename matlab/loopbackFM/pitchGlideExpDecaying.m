% input parameters
fs = 44100;
dur = 0.5;
alpha = 0.00025;
f0_1 = 100;         % starting frequency, Hz
f0_2 = 40;          % ending frequency, Hz
B = 0.01;           % timbre control

% derived parameters
N = dur*fs;
T = 1/fs;
n = 0:1:N-1;
nT = 0:T:(dur-T);

% decay envelope
A_e = 1;
tau_e = -(N-1)/log(0.001);
env = A_e * exp(-n/tau_e);

%w0 = 2*pi*f0;
%w0Tilde = w0*exp(-alpha * n);

w0 = 2*pi*f0_1;
tau = -(N-1)/log(f0_2/f0_1);
w0Tilde = w0 * exp(-n./tau);


%% loopback FM

wc = w0/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0Tilde/wc).^2);

y = zeros(1, N);
y(1) = 1;
for i=2:N
    y(i) = exp(j*wc*T*(1 + BTilde(i) * real(y(i-1)))) * y(i-1);
end

spectrogram(real(y), hann(256), 128, 1024, fs, 'yaxis')
ylim([0 2])

soundsc(real(y) .* env, fs)

%% stretched APF

% exponential decrease for w0
A = w0;
tau = -dur/log(2*pi*40/w0);
b0 = (sqrt(1 - B^2) - 1)/B;
%ThetaH = integral [exp(-nT/tau dn] <-- that's the equation
ThetaH = (-tau * A) * exp(-nT/tau); 
Y = (b0 + exp(1j.*ThetaH)) ./ (1 + b0.*exp(1j.*ThetaH));

% linear decrease for w0
m = 2 * pi * (-60/dur);
b = 2 * pi * 100;
b0 = (sqrt(1 - B^2) - 1)/B;
%ThetaH = integral [m * nT + b d(nT)] <-- that's the equation (why is T
%included here, but not in the version above?
ThetaH = (m*nT.^2)/2 + b*nT;
Y = (b0 + exp(1j.*ThetaH)) ./ (1 + b0.*exp(1j.*ThetaH));

spectrogram(real(Y), hann(256), 128, 1024, fs, 'yaxis')
ylim([0 2])