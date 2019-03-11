% This script tries to emulate the sound of the nonlinear von Karman plate
% from Bilbao's 2005 paper "Sound Synthesis for Nonlinear Plates." It uses
% frequencies that were found using peak picking from the linear version of
% that same plate.

addpath(genpath('/Users/jenniferhsu/Documents/programming/sinemodel'));
addpath(genpath('../proofOfConcept/'));

%% input parameters 

% general
fs = 44100;
dur = 1;
plotSpectrograms = 0;
writeAudioFiles = 0;

% feedback FM
B = 0.5;    % feedback coefficient
B = 0.7;
g = 0.9999; % pitch glide coefficient

[yLin, ~] = audioread('audioExamples/linear.wav');
[yNL, ~] = audioread('audioExamples/nonlinear.wav');

%% derived parameters
N = fs*dur;
T = 1/fs;
n = 0:1:N-1;
nT = 0:T:(dur-T);

%% modal frequencies from linear plate 
% using Dan Ellis's sinemodel package from
% http://www.ee.columbia.edu/~dpwe/LabROSA/matlab/sinemodel/

Nfft = 1024;
S = specgram(yLin, Nfft);               % SP Toolbox routine (or use ifgram.m below)
[R, M] = extractrax(abs(S), 0.0005);	% find peaks in STFT *magnitude*
disp(['size of R is ',num2str(size(R,1)),' rows x ',num2str(size(R,2)),' cols']);
F = R * fs / Nfft;                      % Convert R from bins to Hz

% plot it out
% specgram(yLin, Nfft, fs)
% colormap(1-gray)                        % black is intense, white is quiet
% hold on
% tt = [1:size(F,2)]*(Nfft/2)/fs;         % default specgram step is NFFT/2 i.e. 128
% plot(tt, F', 'r');                      % the tracks follow the specgram peaks
% 
% dr1 = synthtrax(F, M/64, fs, Nfft, Nfft/2); % Divide M by 64 to factor out window, FFT weighting
% specgram(dr1, Nfft, fs)

Nf = size(F, 1);
fVec = zeros(1, Nf);
for i=1:Nf
    freqs = F(i,~isnan(F(i,:)));
    fVec(i) = freqs(1);
end

fVecEnd = fVec/1.08;    % divisor is by eye from Sonic Visualizer

fcVec = fVec/(sqrt(1 - B^2));
fcVecEnd = fVecEnd/(sqrt(1 - B^2));

%% modal frequencies from nonlinear plate 

Nfft = 1024;
SNL = specgram(yNL, Nfft);               % SP Toolbox routine (or use ifgram.m below)
[RNL, MNL] = extractrax(abs(SNL), 0.0005);	% find peaks in STFT *magnitude*
disp(['size of R2 is ',num2str(size(RNL,1)),' rows x ',num2str(size(R,2)),' cols']);
FNL = RNL * fs / Nfft;                      % Convert R from bins to Hz

% plot it out
% specgram(yNL, Nfft, fs)
% colormap(1-gray)                        % black is intense, white is quiet
% hold on
% tt = [1:size(FNL,2)]*(Nfft/2)/fs;         % default specgram step is NFFT/2 i.e. 128
% plot(tt, FNL', 'r');                      % the tracks follow the specgram peaks

% dr1 = synthtrax(FNL, MNL/64, fs, Nfft, Nfft/2); % Divide M by 64 to factor out window, FFT weighting
% specgram(dr1, Nfft, fs)

NfNL = size(FNL, 1);
fVecNL = zeros(1, NfNL);
for i=1:NfNL
    freqs = FNL(i,~isnan(FNL(i,:)));
    fVecNL(i) = freqs(1);
end

NLInds = zeros(1, Nf);
for i=1:Nf
    [~, NLInds(i)] = min(abs(fVec(i) - fVecNL));
end

fVecNL(NLInds)



%% set up decay envelopes (linear)
% grab decay envelopes the sinusoidal analysis

nInds = find(isnan(M));
M(nInds) = 0;
M = M/max(max(M));

env = zeros(Nf, N);

for i=1:Nf
    mEnv = M(i,:);
    
    % get start/end values and indices
    startVal = mEnv(1);
    startInd = 0;
    
    [m, ind] = min(mEnv);
   
    endVal = max(m, .001);
    endInd = ind*(Nfft/2);
    
    % solve the exponentially decaying envelope problem
    A = startVal;
    tau = -endInd/log(endVal/startVal);
    
    env(i,:) = A * exp(-n/tau);
    
end

%% set up decay envelopes (nonlinear, full nonlinear modes)
% grab decay envelopes the sinusoidal analysis

MNL2 = MNL;
nInds = find(isnan(MNL2));
MNL2(nInds) = 0;
MNL2 = MNL2/max(max(MNL2));

envNL1 = zeros(NfNL, N);
xx = linspace(0, size(MNL2, 2)-1, size(MNL2, 2));

for i=1:NfNL
    mEnv = MNL2(i,:);
    p = polyfit(xx,mEnv,3);
    f = polyval(p,xx);
    
    envNL1(i,:) = resample(f, 44100, length(mEnv));
    %envNL1(i,:) = envelope(envNL1(i,:), 4096, 'rms');
    %envNL1(i,:) = envelope(envNL1(i,:), 4096);

end

%% set up decay envelopes (nonlinear, only modes that match linear)
% grab decay envelopes the sinusoidal analysis

MNL2 = MNL(NLInds,:);
nInds = find(isnan(MNL2));
MNL2(nInds) = 0;
MNL2 = MNL2/max(max(MNL2));

envNL = zeros(Nf, N);

for i=1:Nf
    mEnv = MNL2(i,:);
    
    envNL(i,:) = resample(mEnv, 44100, length(mEnv));
    envNL(i,:) = envelope(envNL(i,:), 4096, 'rms');
    
    % get start/end values and indices
%     startVal = max(mEnv);
%     startInd = 0;
%     
%     [m, ind] = min(mEnv);
%    
%     endVal = max(m, .001);
%     endInd = ind*(Nfft/2);
%     
%     % solve the exponentially decaying envelope problem
%     A = startVal;
%     tau = -endInd/log(endVal/startVal);
%     
%     envNL(i,:) = A * exp(-n/tau);
    
end

%% second type of envelope

aStart = zeros(1, Nf);
aStart(1) = 1;              % highest amplitude for lowest modal freq
aStart(Nf) = 0.5;          % lowest amplitude for highest modal freq
tau = -(Nf-1) / log(aStart(Nf));
a0 = 1/exp(-1/tau);
aStart(2:Nf-1) = a0*exp(-(2:Nf-1)/tau);     % exponentially decreasing 
                                            % starting amplitudes

e = g.^(linspace(0, N, N));   % exponential decay

env2 = zeros(Nf, N);
for i=1:Nf
    env2(i,:) = aStart(i) * e;
end


%% modal synthesis with exponential decay using the LINEAR modes
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


%% modal synthesis with exponential decay using the NONLINEAR modes
yMSNL = zeros(1, N);
yMSNLMat = zeros(NfNL, N);
nVec = 0:T:(dur-T);

for i=1:NfNL
    f = fVecNL(i);
    yMSNLMat(i,:) = exp(1j*2*pi*f*nVec) .* envNL1(i,:);
    yMSNL = yMSNL + (exp(1j*2*pi*f*nVec) .* envNL1(i,:));
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
% w0Vec = 2 * pi * fVec';
% tauVec = -(N-1)./log(fVecEnd./fVec)';
% w0VecTilde = w0Vec .* exp(-n./tauVec(:));

A = 1;
tau = -(N-1)./log(0.001);
x = A*exp(-n/tau);
mult = fVec - fVecEnd;
f0VecTilde = zeros(Nf, N);
for i=1:Nf
    f0VecTilde(i,:) = fVecEnd(i) + mult(i) * x;
end
w0VecTilde = 2*pi*f0VecTilde;

wcVec = w0VecTilde(:,1)/sqrt(1 - B^2);
BTilde = sqrt(1 - (w0VecTilde./wcVec).^2);

% envelope from linear model
yLBFM2 = zeros(1, N);
yLBFMMat2 = zeros(Nf, N);
yLBFMMat2(:,1) = 1;
for f=1:Nf
    for i=2:N
        yLBFMMat2(f,i) = exp(j*wcVec(f)*T*(1 + BTilde(f,i) * real(yLBFMMat2(f,i-1)))) * yLBFMMat2(f,i-1);
    end
    yLBFM2 = yLBFM2 + yLBFMMat2(f,:) .* env(f,:);
end

% envelope from nonlinear model
yLBFM2NL = zeros(1, N);
yLBFM2NLMat = zeros(Nf, N);
yLBFM2NLMat(:,1) = 1;
for f=1:Nf
    for i=2:N
        yLBFM2NLMat(f,i) = exp(j*wcVec(f)*T*(1 + BTilde(f,i) * real(yLBFM2NLMat(f,i-1)))) * yLBFM2NLMat(f,i-1);
    end
    yLBFM2NL = yLBFM2NL + yLBFM2NLMat(f,:) .* envNL(f,:);
end

% my own decaying enveloeps
yLBFM22 = zeros(1, N);
yLBFMMat22 = zeros(Nf, N);
yLBFMMat22(:,1) = 1;
for f=1:Nf
    for i=2:N
        yLBFMMat22(f,i) = exp(j*wcVec(f)*T*(1 + BTilde(f,i) * real(yLBFMMat22(f,i-1)))) * yLBFMMat22(f,i-1);
    end
    yLBFM22 = yLBFM22 + yLBFMMat22(f,:) .* env2(f,:);
end

%% Loopback FM using nonlinear modes and envelopes

A = 1;
tau = -(N-1)./log(0.001);
x = A*exp(-n/tau);
fVecEndNL = fVecNL/1.08;
mult = fVecNL - fVecEndNL;
f0VecNLTilde = zeros(NfNL, N);
for i=1:NfNL
    f0VecNLTilde(i,:) = fVecEndNL(i) + mult(i) * x;
end
w0VecNLTilde = 2*pi*f0VecNLTilde;

wcVecNL = w0VecNLTilde(:,1)/sqrt(1 - B^2);
BTildeNL = sqrt(1 - (w0VecNLTilde./wcVecNL).^2);


% envelope from nonlinear model
yLBFM3 = zeros(1, N);
yLBFM3Mat = zeros(NfNL, N);
yLBFM3Mat(:,1) = 1;
for f=1:NfNL
    for i=2:N
        yLBFM3Mat(f,i) = exp(j*wcVecNL(f)*T*(1 + BTildeNL(f,i) * real(yLBFM3Mat(f,i-1)))) * yLBFM3Mat(f,i-1);
    end
    yLBFM3 = yLBFM3 + yLBFM3Mat(f,:) .* envNL1(f,:);
end


%% Time-varying allpass filter tests
TVAPFParams.M = 100;
TVAPFParams.f_m = 500;
TVAPFParams.f_b = 200;

[yTVAPFMemb2, yTVAPFMembMat2, TVAPFParams2] = TVAPFSynthesis(fVec, env, TVAPFParams, 1, fs);

% exponential pitch glide
yTVAPFLB2 = zeros(1, N);
for i=1:Nf
    f_pi = fVec(i);
    %y1 = TVAPF2(yLBFMMat2(i,:), f_pi, TVAPFParams.f_b, TVAPFParams.M, TVAPFParams.f_m, fs);
    %y1 = TVAPF2(yLBFMMat2(i,:), f_pi, TVAPFParams2.f_b(i), TVAPFParams2.M(i), TVAPFParams2.f_m(i), fs);
    %y1 = TVAPF2(yLBFM2NLMat(i,:), f_pi, TVAPFParams.f_b, TVAPFParams.M, TVAPFParams.f_m, fs);
    y1 = TVAPF2(yLBFM2NLMat(i,:), f_pi, TVAPFParams2.f_b(i), TVAPFParams2.M(i), TVAPFParams2.f_m(i), fs);
    yTVAPFLB2 = yTVAPFLB2 + y1 .* env(i,:);
end



%% adding filtered noise

% filter parameters from:
% https://www.cim.mcgill.ca/~clark/nordmodularbook/nm_percussion.html
[B,A] = butter(3, [523 7000]/(fs/2));

wn = 2*rand(1, N) - 1;
fwn = filter(B, A, wn);

yMSN = yMS + fwn .* env(1,:);
yLBFMN1 = yLBFM1 + fwn .* env(1,:);
yLBFMN2 = yLBFM2 + fwn .* x(1,:);

%% maybe do convolutional synthesis instead of filtered noise?
% the pitch glide is pretty on point though

% YES, THEY SOUND PRETTY GOOD.