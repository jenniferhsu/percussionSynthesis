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


%% modal frequencies from nonlinear plate 

Nfft = 1024;
SNL = specgram(yNL, Nfft);               % SP Toolbox routine (or use ifgram.m below)
[RNL, MNL] = extractrax(abs(SNL), 0.0005);	% find peaks in STFT *magnitude*
disp(['size of R2 is ',num2str(size(RNL,1)),' rows x ',num2str(size(RNL,2)),' cols']);
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
fVecNLLast = zeros(1, NfNL);
for i=1:NfNL
    freqs = FNL(i,~isnan(FNL(i,:)));
    fVecNL(i) = freqs(1);
    fVecNLLast(i) = freqs(end);
end

[fVecNL, sortInds] = sort(fVecNL, 'ascend');
fVecNLLast = fVecNLLast(sortInds);
MNLSort = MNL(sortInds, :);
FNLSort = FNL(sortInds, :);


%% set up decay envelopes (nonlinear, full nonlinear modes)
% grab decay envelopes the sinusoidal analysis

MNL2 = MNLSort;
nInds = find(isnan(MNL2));
MNL2(nInds) = 0;
MNL2 = MNL2/max(max(MNL2));

envNL1 = zeros(NfNL, N);
xx = linspace(0, size(MNL2, 2)-1, size(MNL2, 2));

tau = -(N-1)/log(0.1);
ee = 1 * exp(-n/tau);

for i=1:NfNL
    mEnv = MNL2(i,:);

    p = polyfit(xx,mEnv,3);
    f = polyval(p,xx);
    
    envNL1(i,:) = resample(f, 44100, length(mEnv));
%     if i > 1
%         envNL1(i,:) = envNL1(i,:) .* ee; % need to get it to decrease
%     end
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
%for f=1:NfNL
for f=1:(NfNL-60)
    for i=2:N
        yLBFM3Mat(f,i) = exp(j*wcVecNL(f)*T*(1 + BTildeNL(f,i) * real(yLBFM3Mat(f,i-1)))) * yLBFM3Mat(f,i-1);
    end
    yLBFM3 = yLBFM3 + yLBFM3Mat(f,:) .* envNL1(f,:);
end

figure
subplot(211)
spectrogram(yNL, hann(256), 128, 1024, fs, 'yaxis');
subplot(212)
spectrogram(real(yLBFM3), hann(256), 128, 1024, fs, 'yaxis');


%% None of the linear modal frequencies are in the nonlinear modal 
% frequencies list. Maybe in the linear case, wc = w0 and B = 0 and in the
% nonlinear case, wc is fVec (linear), w0 is w0Vec, and B is BVec.

wcVec = 2*pi*fVec;
w0Vec = zeros(1, Nf);
w0VecEnd = zeros(1, Nf);
NLInds = zeros(1, Nf);
for i=1:Nf
    inds = find(fVec(i) - fVecNL > 0);
    [m, ind] = min(fVec(i) - fVecNL(inds));
    
    NLInds(i) = (inds(ind));
    w0Vec(i) = 2*pi*fVecNL(NLInds(i));
    w0VecEnd(i) = 2*pi*fVecNLLast(NLInds(i));
    
    %fprintf('%f, %f, %f\n', m, fVec(i), fVecNL(ind));
end

% This could be what is happening to the the modal frequencies
% when we go from linear to nonlinear. The B starting value is not 1
% (maybe)
BVecStart = sqrt(1 - (w0Vec./wcVec).^2); 

%% Can we synthesize this sce-nare-nare?

% envelope from nonlinear model, no pitch glide
yLBFM4 = zeros(1, N);
yLBFM4Mat = zeros(Nf, N);
yLBFM4Mat(:,1) = 1;
for f=1:Nf
    for i=2:N
        %yLBFM4Mat(f,i) = exp(j*wcVecNL(f)*T*(1 + BVecStart(f) * real(yLBFM4Mat(f,i-1)))) * yLBFM4Mat(f,i-1);
        yLBFM4Mat(f,i) = exp(j*wcVec(f)*T*(1 + BVecStart(f) * real(yLBFM4Mat(f,i-1)))) * yLBFM4Mat(f,i-1);
    end
    yLBFM4 = yLBFM4 + yLBFM4Mat(f,:) .* envNL1(NLInds(f),:);
end

figure
subplot(211)
spectrogram(yNL, hann(256), 128, 1024, fs, 'yaxis');
subplot(212)
spectrogram(real(yLBFM4), hann(256), 128, 1024, fs, 'yaxis');


% so our next step is to create a pitch glide using these
% let's make it similar to that pitch glide according to w0Vec and w0VecEnd

f0V = w0Vec/(2*pi);
f0VEnd = w0VecEnd/(2*pi);

A = 1;
tau = -(N*T)./log(0.001);
x = A*exp(-nT/tau);

% the pitch glides
f0Exp = (f0VEnd(:) - f0V(:)) * x + f0V(:);

% synthesis with LBFM
w0VTilde = 2*pi*f0Exp;
wcVTilde = zeros(Nf, N);
for f=1:Nf
    wcVTilde(f,:) = w0VTilde(f,1)/sqrt(1 - BVecStart(f).^2);
end
BTilde = sqrt(1 - (w0VTilde./wcVTilde).^2);
yLBFM5 = zeros(1, N);
yLBFM5Mat = zeros(NfNL, N);
yLBFM5Mat(:,1) = 1;
%for f=1:NfNL
for f=1:Nf
    for i=2:N
        yLBFM5Mat(f,i) = exp(j*wcVTilde(f)*T*(1 + BTilde(f,i) * real(yLBFM5Mat(f,i-1)))) * yLBFM5Mat(f,i-1);
    end
    yLBFM5 = yLBFM5 + yLBFM5Mat(f,:) .* envNL1(NLInds(f),:);
end



%% Let's see if we can find any modulation frequencies that are common in 
% the nonlinear frequencies list
diffMat = zeros(NfNL, NfNL);
thresh = 1;

for i=1:NfNL
    diffMat(i,:) = abs(fVecNL(i) - fVecNL);
end

% this will create a list of frequencies that show up in the differences
% matrix
modF = [];
modIndsI = [];
modIndsK = [];
for ind=1:NfNL
    for i=1:NfNL
        for k=i+1:NfNL
            if(abs(diffMat(i,k) - fVecNL(ind)) <= thresh)
                modF = [modF fVecNL(ind)];
                modIndsI = [modIndsI i];
                modIndsK = [modIndsK k];
            end
        end
    end
end

modFI = zeros(size(modF));
modFK = zeros(size(modF));
for i=1:length(modF)
    modFI(i) = fVecNL(modIndsI(i));
    modFK(i) = fVecNL(modIndsK(i));
end

% idea, check if we can find these kinds of modulations using the linear
% modes from the Bilbao plate

% so now:
% abs(modFI' - modFK') is equivalent to modF'
% and modF holds frequencies that are in fVecNL

[sortModF, sInds] = sort(modF, 'ascend');
sortModFI = modFI(sInds);
sortModFK = modFK(sInds);

% if i make proper modulations, i can take out the higher, larger
% frequencies and replace with modulating frequencies

% sortModF' + sortModFI' OR sortModF' + sortModFK'
% will create the higher frequency that I can take out
% so things in sortModF MIGHT be in the list of linear modal frequencies


