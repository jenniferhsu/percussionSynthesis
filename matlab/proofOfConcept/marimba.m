% This file is an attempt to make the sound of a marimba (or a bar with
% two free ends)

% (page 47 of the Science of Percussion Instruments)

%% input parameters 

% general
fs = 44100;
dur = 2;
plotSpectrograms = 0;
writeAudioFiles = 1;
outDir = 'audioExamples/marimba/';

if ~exist(outDir)
    mkdir(outDir)
end

% feedback FM
B = 0.3;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

% time-varying APF
TVAPFParams.M = fs/40;
TVAPFParams.f_m = 100;
TVAPFParams.f_b = fs/16;

%% derived parameters
N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient

%% ideal bar modal frequencies
% from Science of Percussion Instruments page 47

f1 = 440;
Nf = 7;
fVecBar = zeros(1, Nf);
for i=1:Nf
    if i==1
        fVecBar(i) = 3.011^2;
    else
        fVecBar(i) = (2*i+1)^2;
    end
end
fVecBar = fVecBar/fVecBar(1);
fVecBar = f1 * fVecBar;


% center frequencies for feedback FM if sounding frequencies = membrane
% modal freqeuncies
fcVecBar = fVecBar./(sqrt(1-B^2));

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

e = g.^(linspace(0, N, N));   % exponential decay

env = zeros(Nf, N);
for i=1:Nf
    env(i,:) = aStart(i) * e;
end

%% modal synthesis with exponential decay using the modes
yMS = zeros(1, N);
nVec = 0:T:(dur-T);

for i=1:Nf
    f = fVecBar(i);
    yMS = yMS + (exp(1j*2*pi*f*nVec) .* env(i,:));
end

% aStart(1:Nf) = ones(1, Nf);
% yMS2 = zeros(1, N);
% for i=1:Nf
%     f = fVecBar(i);
%     yMS2 = yMS2 + (exp(1j*2*pi*f*nVec) .* (aStart(i)*e));
% end

if plotSpectrograms == 1
    figure
    spectrogram(real(yMS), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
end


%% FEEDBACK FM

% EXAMPLE 1: feedback FM center frequencies = membrane modal frequencies
% (YES)
[yFBFMMemb1, yFBFMMembMat1] = feedbackFMSynthesis(fcVecBar, B, env, fs);

% EXAMPLE 2: feedback FM with pitch glide with center frequencies = membrane
% modal frequencies 
BVec = linspace(0.3, 0.4, N);
[yFBFMMemb2, yFBFMMembMat2] = feedbackFMSynthesis(fcVecBar, BVec, env, fs);

% EXAMPLE 3:
% for plotting purposes
BVec2 = linspace(0.5, 0.7, N);
fcVecBar2 = fVecBar./(sqrt(1-BVec2(1)^2));
[yFBFMMemb3, yFBFMMembMat3] = feedbackFMSynthesis(fcVecBar2, BVec2, env, fs);

% this is the pitch glide function
%w0 = (2*pi*fVecBar(1))*sqrt(1 - linspace(0.9, 0.91, N).^2);
%plot(w0/(2*pi));

%% STRETCHED APF

b0 = (sqrt(1-B^2) - 1)/B;
fVecBar_w0 = fcVecBar*sqrt(1 - B^2);

% EXAMPLE 1: w0 = modal frequencies 
[ySAPFMemb1, ySAPFMembMat1] = stretchedAPFSynthesis(fVecBar, b0, env, fs, [], 'none');

% EXAMPLE 2: w0 = modal frequencies with a pitch glide
% this pitch glide sounds nice
[ySAPFMemb2, ySAPFMembMat2] = stretchedAPFSynthesis(fVecBar, b0, env, fs, fVecBar - 20*ones(1, Nf), 'linear');

% EXAMPLE 3: pitch glide from FBFM Example 2 - this is really well done
fVecBar_w0End = fVecBar .* sqrt(1 - BVec(end).^2);
[ySAPFMemb3, ySAPFMembMat3] = stretchedAPFSynthesis(fVecBar, b0, env, fs, fVecBar_w0End, 'linear');


%% writing to file
if writeAudioFiles == 1
    
    audiowrite([outDir 'yMS.wav'], scaleForSavingAudio(real(yMS)), fs);
    
    audiowrite([outDir 'yFBFMMemb1.wav'], scaleForSavingAudio(real(yFBFMMemb1)), fs);
    audiowrite([outDir 'yFBFMMemb2.wav'], scaleForSavingAudio(real(yFBFMMemb2)), fs);
    audiowrite([outDir 'yFBFMMemb3.wav'], scaleForSavingAudio(real(yFBFMMemb3)), fs);
    
    audiowrite([outDir 'ySAPFMemb1.wav'], scaleForSavingAudio(real(ySAPFMemb1)), fs);
    audiowrite([outDir 'ySAPFMemb2.wav'], scaleForSavingAudio(real(ySAPFMemb2)), fs);
    audiowrite([outDir 'ySAPFMemb3.wav'], scaleForSavingAudio(real(ySAPFMemb3)), fs);
end