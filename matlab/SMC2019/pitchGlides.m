% Pitch glide with the stretched allpass filter and feedback FM
% formulations for exponential, linear, and logarithmic B(n)
fs = 44100;
dur = 1;

T = 1/fs;
N = fs * dur;
nT = 0:T:(dur-T);

wc = 2*pi*698; 
g = 0.9999;

C = 0;                      % constant of integration

%% exponential B(n)
% B = g^n
%g = 0.99999; % plot
BExp = g.^(0:N-1)';
%BExp = 1/(g.^(0:N-1)');

% this one is a good one to look at
BExp = sqrt(sqrt(g.^(0:N-1)));

% the expected pitch trajectory
w0Exp = wc * sqrt(1-BExp.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stretched Allpass Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = sqrt(1 - BExp.^2);
ThetaExp = (wc*T/log(g)) .* (u - atanh(u) + C);

b0Exp = (sqrt(1 - BExp.^2) - 1)./BExp;

H_apf_exp = (b0Exp + exp(1j.*ThetaExp)) ./ (1 + b0Exp.*exp(1j.*ThetaExp));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback FM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_fm_exp = zeros(1, N);
h_fm_exp(1) = 1;

for n=2:N
    h_fm_exp(n) = (exp(1j*wc*T*(1 + BExp(n) * real(h_fm_exp(n-1))))) * h_fm_exp(n-1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback PM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_pm_exp = zeros(1, N);
h_pm_exp(1) = 1;

for n=1:N-1
    ind = n+1;
    h_pm_exp(ind) = exp(1j*wc*T*(n + real(sum(BExp(1:n)' .* h_pm_exp(1:ind-1)))));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(311)
spectrogram(H_apf_exp, hann(1024), 512, 1024, fs, 'yaxis')
hold on
plot(linspace(0, 1000, N), w0Exp/(2*pi*1000), 'r')
ylim([0 2])
title('Allpass Filter Pitch Glide with exponential B(n) and w0(n)');
legend('w_0(n)')

subplot(312)
spectrogram(h_fm_exp, hann(1024), 512, 1024, fs, 'yaxis')
hold on
plot(linspace(0, 1000, N), w0Exp/(2*pi*1000), 'r')
ylim([0 2])
title('Feedback FM Pitch Glide with exponential B(n) and w0(n)');
legend('w_0(n)')

subplot(313)
spectrogram(h_pm_exp, hann(1024), 512, 1024, fs, 'yaxis')
hold on
plot(linspace(0, 1000, N), w0Exp/(2*pi*1000), 'r') 
ylim([0 2])
title('Feedback PM Pitch Glide with exponential B(n) and w0(n)');
legend('w_0(n)')


%% linear B(n)
% B(n) = n;

BLin = nT';
%BLin = sin(2*pi*100*nT');

% the expected pitch trajectory
w0Lin = wc * sqrt(1-BLin.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stretched Allpass Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ThetaLin = wc*0.5*nT'.*sqrt(1 - BLin.^2) + wc*0.5*asin(BLin) + C;

b0Lin = (sqrt(1 - BLin.^2) - 1)./BLin;

H_apf_lin = (b0Lin + exp(1j.*ThetaLin)) ./ (1 + b0Lin.*exp(1j.*ThetaLin));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback FM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_fm_lin = zeros(1, N);
h_fm_lin(1) = 1;

for n=2:N
    h_fm_lin(n) = (exp(1j*wc*T*(1 + BLin(n) * real(h_fm_lin(n-1))))) * h_fm_lin(n-1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback PM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_pm_lin = zeros(1, N);
h_pm_lin(1) = 1;

for n=1:N-1
    ind = n+1;
    h_pm_lin(ind) = exp(1j*wc*T*(n + real(sum(BLin(1:n)' .* h_pm_lin(1:ind-1)))));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

subplot(311)
spectrogram(H_apf_lin, hann(1024), 512, 1024, fs, 'yaxis')
hold on
plot(linspace(0, 1000, N), w0Lin/(2*pi*1000), 'r')
ylim([0 2])
title('Allpass Filter Pitch Glide with linear B(n) and w0(n)');
legend('w_0(n)')

subplot(312)
spectrogram(h_fm_lin, hann(1024), 512, 1024, fs, 'yaxis')
hold on
plot(linspace(0, 1000, N), w0Lin/(2*pi*1000), 'r')
ylim([0 2])
title('Feedback FM Pitch Glide with linear B(n) and w0(n)');
legend('w_0(n)')

subplot(313)
spectrogram(h_pm_lin, hann(1024), 512, 1024, fs, 'yaxis')
hold on
plot(linspace(0, 1000, N), w0Lin/(2*pi*1000), 'r')  % the expected linear sweep
ylim([0 2])
title('Feedback PM Pitch Glide with linear B(n) and w0(n)');
legend('w_0(n)')


%% linear w0(n)

m = 780/dur;
b = 220;
f0 = m*nT + b;
w0 = 2*pi*f0;

% linear w0(n) means that B is of this form:
B = sqrt((1 - (m*nT + b).^2)/wc);

Theta = (m/2)*nT.^2 + b*nT + C;

b0 = (sqrt(1 - B.^2) - 1)./B;

H_apf = (b0 + exp(1j.*Theta)) ./ (1 + b0.*exp(1j.*Theta));

figure
spectrogram(H_apf, hann(1024), 512, 1024, fs, 'yaxis')
hold on
plot(linspace(0, 1000, N), w0/(2*pi*1000), 'r')
ylim([0 2])
title('Allpass Filter Pitch Glide with linear B(n) and w0(n)');
legend('w_0(n)')

% this is kind of how it should sound
x1 = sin(2*pi*cumsum(w0/(2*pi))*T);

