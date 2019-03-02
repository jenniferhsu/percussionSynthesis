% TVAPF2 parameter section plotting

addpath(genpath('../proofOfConcept/'));

fs = 44100;
dur = 1;

N = fs*dur;
T = 1/fs;
nT = 0:T:(dur-T);

% FFT parameters
Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

%% Sinusoidal input

f = 5500;
x = sin(2*pi*f*nT);

% TVAPF parameters
f_pi = f;   %   f_pi: the frequency at which the phase is -pi
f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 1000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

y = TVAPF2(x, f_pi, f_b, M, f_m, fs);
Y = fft(y, Nfft);
YPos = Y(1:Nfft/2+1);
YPosAbsdB = 20*log10(abs(YPos)/max(abs(YPos)));

figure
plot(faxis, YPosAbsdB, 'linewidth', 2)
hold on
%plot([f f], [-60 0], 'r--');
% for i=1:5
%     fout1 = f+(i*f_m);
%     % aliased sidebands
%     if fout1 > fs/2
%         fout1 = (fs/2) - (fout - fs/2);
%     end
%     plot([fout1 fout1], [-60 0], 'g--');
%     
%     fout2 = f-(i*f_m);
%     % aliased sidebands
%     if fout2 < 0
%         fout2 = -1*fout2;
%     end
%     plot([fout2 fout2], [-60 0], 'm--');
% end
%xlim([faxis(1) faxis(end)])
ylim([-60 0])
grid on
fig = gcf
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2.5];
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Sinusoid filtered by a time-varying allpass mag spectrum');
print('figures/timeVaryingAPFSineInput', '-depsc', '-r0')


% figure
% spectrogram(y, hann(256), 128, 1024, fs, 'yaxis');

%% Sinusoidal Input (Aliasing)
f = 300;
x = sin(2*pi*f*nT);

% TVAPF parameters
f_pi = f;   %   f_pi: the frequency at which the phase is -pi
f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 1000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

yAlias = TVAPF2(x, f_pi, f_b, M, f_m, fs);
YAlias = fft(yAlias, Nfft);
YAliasPos = YAlias(1:Nfft/2+1);
YAliasPosAbsdB = 20*log10(abs(YAliasPos)/max(abs(YAliasPos)));

figure
plot(faxis, YAliasPosAbsdB, 'linewidth', 2)
hold on
%plot([f f], [-60 0], 'r--');
% for i=1:5
%     fout1 = f+(i*f_m);
%     % aliased sidebands
%     if fout1 > fs/2
%         fout1 = (fs/2) - (fout - fs/2);
%     end
%     plot([fout1 fout1], [-60 0], 'g--');
%     
%     fout2 = f-(i*f_m);
%     % aliased sidebands
%     if fout2 < 0
%         fout2 = -1*fout2;
%     end
%     plot([fout2 fout2], [-60 0], 'm--');
% end
%xlim([faxis(1) faxis(end)])
ylim([-60 0])
grid on
fig = gcf
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2.5];
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Sinusoid filtered by a time-varying allpass mag spectrum (with aliasing)');
print('figures/timeVaryingAPFSineInputAlias', '-depsc', '-r0')

%% Static Loopback FM input
f = 700;
w0 = 2*pi*f;
B = 0.9;
wc = w0/sqrt(1 - B^2);

% loopback FM
h = zeros(1, N);
h(1) = 1;
for n=2:N
    h(n) = exp(1j*wc*T*(1 + B*real(h(n-1)))) * h(n-1);
end

% fft analysis
H = fft(h, Nfft);
HPos = H(1:Nfft/2+1);
HPosAbsdB = 20*log10(abs(HPos)/max(abs(HPos)));

% TVAPF
f_pi = f;   %   f_pi: the frequency at which the phase is -pi
f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 1000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

yLB = TVAPF2(real(h), f_pi, f_b, M, f_m, fs);

% fft analysis
YLB = fft(yLB, Nfft);
YLBPos = YLB(1:Nfft/2+1);
YLBPosAbsdB = 20*log10(abs(YLBPos)/max(abs(YLBPos)));

% TV APF predicted peaks for f at 700
K = 10;
aliasedPeakFreqs = zeros(K, 32);
for k=1:K
    for i=1:2:16
        aliasedPeakFreqs(k,i) = k*f + i*f_m;
        if aliasedPeakFreqs(k,i+1) > (fs/2)
            aliasedPeakFreqs(k,i+1) = (fs/2) - (aliasedPeakFreqs(k,i+1) - (fs/2));
        end
        
        aliasedPeakFreqs(k,i+1) = k*f - i*f_m;
        if aliasedPeakFreqs(k,i+1) < 0
            aliasedPeakFreqs(k,i+1) = -1 * aliasedPeakFreqs(k,i+1);
        end
    end
end

figure
subplot(211)
plot(faxis, HPosAbsdB, 'linewidth', 2)
hold on
%plot(faxis, HPosAbsdB)
% loopback FM predicted peaks
for i=1:16
    plot(i*[f f], [-60 0], 'r--');
end
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Static loopback FM');
ylim([-60 0])
grid on
% TV APF predicted peaks for f at 700
% for k=1:K
%     for i=1:length(aliasedPeakFreqs(k,:))
%         plot([aliasedPeakFreqs(k,i) aliasedPeakFreqs(k,i)], [-60 0], 'g--');
%     end
% end
subplot(212)
plot(faxis, YLBPosAbsdB, 'linewidth', 2)
hold on
% loopback FM predicted peaks
for i=1:16
    plot(i*[f f], [-60 0], 'r--');
end
ylim([-60 0])
grid on
% fig = gcf
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 6 2.5];
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Static loopback FM filtered by a time-varying allpass mag spectrum');
print('figures/timeVaryingAPFStaticLoopbackFM', '-depsc', '-r0')


%% Are these the same results as if I TV APF'd a sinusoid at 700Hz
% then a sinusoid at 1400Hz
% then a sinusoid at 2100Hz
% and then added them together??

nT = (0:N-1)*T;
amps = [10960 7503 5165 2779 2070 1223 756.7 515.6 280.4 200.1];
amps = amps/max(amps);

f_pi = f;   %   f_pi: the frequency at which the phase is -pi
f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 1000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

xSineMat = zeros(10, N);
ySineAPF = zeros(10, N);
YSineAPFPos = zeros(10, Nfft/2+1);

ySineAPFOut = zeros(1, N);

for k=1:K
    % make sinusoids
    xSineMat(k,:) = sin(2*pi*(k*f)*nT) * amps(k);
    
    % TVAPF the sinusoids
    ySineAPF(k,:) = TVAPF2(xSineMat(k,:), f_pi, f_b, M, f_m, fs);
    
    % add them together
    ySineAPFOut = ySineAPFOut + ySineAPF(k,:);
    
    % fft processing
    Y = fft(ySineAPF(k,:), Nfft);
    YSineAPFPos(k,:) = Y(1:Nfft/2+1);
    
end

YSineAPFOut = fft(ySineAPFOut, Nfft);
YSineAPFOutPos = YSineAPFOut(1:Nfft/2+1);
YSineAPFOutPosdB = 20*log10(abs(YSineAPFOutPos)/max(abs(YSineAPFOutPos)));

figure
plot(faxis, YSineAPFOutPosdB, 'linewidth', 2)
hold on
%plot(faxis, HPosAbsdB)
% loopback FM predicted peaks
for i=1:16
    plot(i*[f f], [-60 0], 'r--');
end
% TV APF predicted peaks for f at 700
% for k=1:K
%     for i=1:length(aliasedPeakFreqs(k,:))
%         plot([aliasedPeakFreqs(k,i) aliasedPeakFreqs(k,i)], [-60 0], 'g--');
%     end
% end
ylim([-60 0])
grid on
fig = gcf
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2.5];
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Sinusoids tuned to loopback FM filtered by a time-varying allpass mag spectrum');




%% Static stretched allpass input -- looks the same as the feedback FM one
nT = (0:N-1)*T;
b0 = (sqrt(1 - B^2) - 1)/B;
H = (b0 + exp(1j*w0*nT)) ./ (1 + b0.*exp(1j*w0*nT));

% TVAPF parameters
f_pi = f;   %   f_pi: the frequency at which the phase is -pi
f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 1000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

ySAPF = TVAPF2(real(H), f_pi, f_b, M, f_m, fs);
YSAPF = fft(ySAPF, Nfft);
YSAPFPos = YSAPF(1:Nfft/2+1);
YSAPFPosAbsdB = 20*log10(abs(YSAPFPos)/max(abs(YSAPFPos)));

figure
plot(faxis, YSAPFPosAbsdB, 'linewidth', 2)
hold on
ylim([-60 0])
grid on
fig = gcf
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2.5];
%print('figures/timeVaryingAPFStaticSAPF', '-depsc', '-r0')


%% Loopback FM with time-varying pitch and timbre input
f = 700;
w0 = 2*pi*f;
B = 0.9;
wc = w0/sqrt(1 - B^2);

g = 0.9999;
BExp = g.^(0:N-1)';

htv = zeros(1, N);
htv(1) = 1;

for n=2:N
    htv(n) = (exp(1j*wc*T*(1 + BExp(n) * real(htv(n-1))))) * htv(n-1);
end

% fft analysis
HTV = fft(real(htv), Nfft);
HTVPos = HTV(1:Nfft/2+1);
HTVPosAbsdB = 20*log10(abs(HTVPos)/max(abs(HTVPos)));

% TVAPF
f_pi = f;   %   f_pi: the frequency at which the phase is -pi
f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 1000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

yLBTV = TVAPF2(real(htv), f_pi, f_b, M, f_m, fs);

% fft analysis
YLBTV = fft(yLBTV, Nfft);
YLBTVPos = YLBTV(1:Nfft/2+1);
YLBTVPosAbsdB = 20*log10(abs(YLBTVPos)/max(abs(YLBTVPos)));

% TV APF predicted peaks for f at 700
K = 10;
aliasedPeakFreqs = zeros(K, 32);
for k=1:K
    for i=1:2:16
        aliasedPeakFreqs(k,i) = k*f + i*f_m;
        if aliasedPeakFreqs(k,i+1) > (fs/2)
            aliasedPeakFreqs(k,i+1) = (fs/2) - (aliasedPeakFreqs(k,i+1) - (fs/2));
        end
        
        aliasedPeakFreqs(k,i+1) = k*f - i*f_m;
        if aliasedPeakFreqs(k,i+1) < 0
            aliasedPeakFreqs(k,i+1) = -1 * aliasedPeakFreqs(k,i+1);
        end
    end
end

figure
plot(faxis, YLBTVPosAbsdB, 'linewidth', 2)
hold on
ylim([-60 0])
grid on
fig = gcf
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2.5];
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Static loopback FM filtered by a time-varying allpass mag spectrum');


figure
subplot(211)
spectrogram(htv, hann(256), 128, 1024, fs, 'yaxis');
title('loopback FM spectrogram')
set(gca, 'fontsize', 15)
ylim([0 18])
subplot(212)
spectrogram(yLBTV, hann(256), 128, 1024, fs, 'yaxis');
set(gca, 'fontsize', 15)
title('time-varying allpass filtered loopback FM spectrogram')
set(gca, 'fontsize', 15)
ylim([0 18])
print('figures/timeVaryingAPFTVLoopbackFM', '-depsc', '-r0')