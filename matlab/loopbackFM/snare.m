%% Let's try to make a snare sound using handpicked frequencies

addpath(genpath('../proofOfConcept/'));

%% input parameters 

% general
fs = 44100;
dur = 0.6;
plotSpectrograms = 0;
writeAudioFiles = 0;

% feedback FM
B = 0.5;    % feedback coefficient
B = 0.5;
g = 0.9999; % pitch glide coefficient

%% derived parameters
N = fs*dur;
T = 1/fs;
n = 0:1:N-1;
nT = 0:T:(dur-T);

%% snare frequencies

% this is handpicked using Sonic Visualizer and the first Snare recording
% in Zyskowski recordings
%fVec = [140.625, 234.375, 281.25, 421.875, 703.125, 890.650];
%fVecEnd = fVec/1.7;

% this one kind of works for loopback FM
%fVec = [185, 330, 747, 2090, 3845, 4464, 7278, 10518, 11130, 15423];
%fVecEnd = fVec/1.7;

% these two frequencies are from:
% https://www.cim.mcgill.ca/~clark/nordmodularbook/nm_percussion.html
fVec = [185, 330];
fVecEnd = [140, 280];       % pitch glide end vector (by eye)

Nf = length(fVec);

fcVec = fVec/(sqrt(1 - B^2));
fcVecEnd = fVecEnd/(sqrt(1 - B^2));

%% set up decay envelopes
% want the lowest modal frequencies to have a larger starting amplitude 
% than the highest modal frequencies

aStart = zeros(1, Nf);
aStart(1) = 1;              % highest amplitude for lowest modal freq
aStart(Nf) = 0.01;          % lowest amplitude for highest modal freq

tau = -(Nf-1) / log(aStart(Nf));
a0 = 1/exp(-1/tau);
aStart(2:Nf-1) = a0*exp(-(2:Nf-1)/tau);     % exponentially decreasing 
                                            % starting amplitudes

%e = g.^(linspace(0, N, N));   % exponential decay
e = exp(-nT/(-dur/log(0.001)));

env = zeros(Nf, N);
for i=1:Nf
    env(i,:) = aStart(i) * e;
end

%% modal synthesis with exponential decay using the modes
yMS = zeros(1, N);
yMSMat = zeros(Nf, N);
nVec = 0:T:(dur-T);

for i=1:Nf
    f = fVec(i);
    yMSMat(i,:) = exp(1j*2*pi*f*nVec) .* env(i,:);
    yMS = yMS + (exp(1j*2*pi*f*nVec) .* env(i,:));
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yMS), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
end

%% LOOPBACK FM

% loopback FM sounding frequencies = snare modal frequencies
[yLBFM1, yLBFMMat1] = feedbackFMSynthesis(fcVec, B, env, fs);

% loopback FM with exponential pitch glide
w0Vec = 2 * pi * fVec';
tauVec = -(N-1)./log(fVecEnd./fVec)';

w0VecTilde = w0Vec .* exp(-n./tauVec(:));
wcVec = w0Vec/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0VecTilde./wcVec).^2);

yLBFM2 = zeros(1, N);
yLBFMMat2 = zeros(Nf, N);
yLBFMMat2(:,1) = 1;
for f=1:Nf
    for i=2:N
        yLBFMMat2(f,i) = exp(j*wcVec(f)*T*(1 + BTilde(f,i) * real(yLBFMMat2(f,i-1)))) * yLBFMMat2(f,i-1);
    end
    yLBFM2 = yLBFM2 + yLBFMMat2(f,:) .* env(f,:);
end

% loopback FM with linear pitch glide
mVec = 2 * pi * ((fVecEnd - fVec)/dur);
bVec = 2 * pi * fVec;
w0VecTilde = mVec(:).*nT + bVec(:);
wcVec = w0Vec/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0VecTilde./wcVec).^2);

yLBFM3 = zeros(1, N);
yLBFMMat3 = zeros(Nf, N);
yLBFMMat3(:,1) = 1;
for f=1:Nf
    for i=2:N
        yLBFMMat3(f,i) = exp(j*wcVec(f)*T*(1 + BTilde(f,i) * real(yLBFMMat3(f,i-1)))) * yLBFMMat3(f,i-1);
    end
    yLBFM3 = yLBFM3 + yLBFMMat3(f,:) .* env(f,:);
end

%% Time-varying allpass filter tests
TVAPFParams.M = 1500;
TVAPFParams.f_m = 100;
TVAPFParams.f_b = 2750;

% exponential pitch glide
yTVAPFLB2 = zeros(1, N);
for i=1:Nf
    f_pi = fVec(i);
    y1 = TVAPF2(yLBFMMat2(i,:), f_pi, TVAPFParams.f_b, TVAPFParams.M, TVAPFParams.f_m, fs);
    yTVAPFLB2 = yTVAPFLB2 + y1 .* env(i,:);
end

% linear pitch glide
yTVAPFLB3 = zeros(1, N);
for i=1:Nf
    f_pi = fVec(i);
    y1 = TVAPF2(yLBFMMat3(i,:), f_pi, TVAPFParams.f_b, TVAPFParams.M, TVAPFParams.f_m, fs);
    yTVAPFLB3 = yTVAPFLB3 + y1 .* env(i,:);
end

%% randomly selected time-varying allpass filter parameters

[yTVAPFMS, ~, TVAPFParams1] = TVAPFSynthesis(fVec, env, TVAPFParams, 1, fs);

% exponential pitch glide
yTVAPFLB2 = zeros(1, N);
for i=1:Nf
    f_pi = fVec(i);
    y1 = TVAPF2(yLBFMMat2(i,:), f_pi, TVAPFParams1.f_b(i), TVAPFParams1.M(i), TVAPFParams1.f_m(i), fs);
    yTVAPFLB2 = yTVAPFLB2 + y1 .* env(i,:);
end

% linear pitch glide
yTVAPFLB3 = zeros(1, N);
for i=1:Nf
    f_pi = fVec(i);
    y1 = TVAPF2(yLBFMMat3(i,:), f_pi, TVAPFParams1.f_b(i), TVAPFParams1.M(i), TVAPFParams1.f_m(i), fs);
    yTVAPFLB3 = yTVAPFLB3 + y1 .* env(i,:);
end

%% adding filtered noise

% filter parameters from:
% https://www.cim.mcgill.ca/~clark/nordmodularbook/nm_percussion.html
[B,A] = butter(3, [523 7000]/(fs/2));

wn = 2*rand(1, N) - 1;
fwn = filter(B, A, wn);

yMSN = yMS + fwn .* env(1,:);
yLBFMN1 = yLBFM1 + fwn .* env(1,:);
yLBFMN2 = yLBFM2 + fwn .* env(1,:);
yLBFMN3 = yLBFM3 + fwn .* env(1,:);

% YES, THEY SOUND PRETTY GOOD.


%% plots
figure
spectrogram(real(yMSN), hann(256), 128, 1024, fs, 'yaxis');
ylim([0 10])

figure
subplot(311)
spectrogram(real(yLBFMN1), hann(256), 128, 1024, fs, 'yaxis');
ylim([0 10])
subplot(312)
spectrogram(real(yLBFMN2), hann(256), 128, 1024, fs, 'yaxis');
ylim([0 10])
subplot(313)
spectrogram(real(yLBFMN3), hann(256), 128, 1024, fs, 'yaxis');
ylim([0 10])