fs = 44100;
dur = 1.0;
w0 = 0.88;
%w60 = 0.001;
T60 = 0.65;
n60 = T60*fs;

w60 = w0 / 10^(60/20);

N = fs*dur;
T = 1/fs;
nT = linspace(0, dur, N);

% w = A * exp(-nT/tau)

A = w0;
tau = (-n60 * T)/log(w60/A);

w = A * exp(-nT/tau);


figure
plot(nT, 20*log10(abs(w)/max(abs(w))))
hold on
plot([0 nT(end)], [-60 -60], 'r--')
plot([T60 T60], [-100 0], 'r--')

% 20*log10(A0) - 20*log10(A60) = 60
% 20(log10(A0) - log10(A60)) = 60
% log10(A0) - log10(A60) = 60/20
% log10(A0/A60) = 60/20
% A0/A60 = 10^(60/20)
% A60 = A0/10^(60/20)
