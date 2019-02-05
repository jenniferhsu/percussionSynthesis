% plot excitations (raised cosine and filtered noise burst for SMC 2019
% paper)

fs = 44100;
dur = 1;
N = fs*dur;

% raised cosine
winLength = 1024;
e_rc = zeros(N,1); 
n = winLength/2:winLength-1;
w = 0.5 * (1 - cos((2*pi*n)/(winLength-1)));

e_rc(1:winLength/2) = ones(winLength/2, 1);
e_rc(winLength/2+1:winLength) = w;

% filtered noise burst

durNB = 0.03;
flow = 200;
fhigh = 2000;

% noise burst
sampNB = ceil(durNB*fs);
e_nb = [(2*rand(1, sampNB))-1 zeros(1, N-sampNB)]';

% filter
[B, A] = butter(2, [flow/(fs/2) fhigh/(fs/2)], 'bandpass');
e_nb = filter(B, A, e_nb);

figure

subplot(211)
plot(e_rc, 'linewidth', 2);
xlabel('Time (samples)');
ylabel('Amplitude (linear)');
title('Raised cosine excitation, window length L=1024 samples');
xlim([0 2048]);
ylim([0 1.1]);
grid on
set(gca, 'FontSize', 15);

subplot(212)
plot(e_nb, 'linewidth', 2);
xlabel('Time (samples)');
ylabel('Amplitude (linear)');
title('Filtered noise burst excitation, \tau _d = 0.01 sec, fc_{low} = 200, fc_{high} = 2000');
xlim([0 2048]);
ylim([-0.8 0.8]);
grid on
set(gca, 'FontSize', 15);

saveas(gcf, 'figures/excitations', 'epsc');




