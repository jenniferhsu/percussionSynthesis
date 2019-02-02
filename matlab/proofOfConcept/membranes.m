% This file is an attempt to make lots of different kinds of sounds

%% input parameters 

% general
fs = 44100;
dur = 1;
plotSpectrograms = 0;
writeAudioFiles = 1;

% feedback FM
B = 0.9;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

% time-varying APF
TVAPFParams.M = fs/40;
TVAPFParams.f_m = 100;
TVAPFParams.f_b = fs/16;

%% derived parameters
N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient

%% ideal membrane modal frequencies
% from Science of Percussion Instruments

% KETTLEDRUM (page 8)
%f1 = 150; % fundamental frequency of a kettledrum
%fVecMembrane = f1 * [0.63, 1, 1.34, 1.44, 1.66, 1.83, 1.98, 2.20, 2.26, 2.29, 2.55, 2.61, 2.66, 2.89];
%outDir = 'audioExamples/membrane/kettledrum/';

% SNARE (page 29)
% a lot of the techniques are too pitched for this to sound good. the ones
% that sound nice are the randomized TVAPF and mode 4 for the stretched APF
% and TVAPF
%fVecMembrane = [299, 331, 507, 616, 674, 582, 859];
%outDir = 'audioExamples/membrane/snare/';

% BASS DRUM (page 35)
%fVecMembrane = [44, 104, 76, 82, 120, 160, 198, 240];
%outDir = 'audioExamples/membrane/bass-drum/';

% O-DAIKO (page 43)
fVecMembrane = [202, 227, 333, 383, 468, 492, 544, 621, 695, 739, 865, 905, 1023];
outDir = 'audioExamples/membrane/o-daiko/';

Nf = length(fVecMembrane);

% center frequencies for feedback FM if sounding frequencies = membrane
% modal freqeuncies
fcVecMembrane = fVecMembrane./(sqrt(1-B^2));

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
    f = fVecMembrane(i);
    yMS = yMS + (exp(1j*2*pi*f*nVec) .* env(i,:));
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yMS), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
end


%% FEEDBACK FM

% EXAMPLE 1: feedback FM center frequencies = membrane modal frequencies
[yFBFMMemb1, yFBFMMembMat1] = feedbackFMSynthesis(fVecMembrane, B, env, fs);

% EXAMPLE 2: feedback FM sounding frequencies = membrane modal frequencies
[yFBFMMemb2, yFBFMMembMat2] = feedbackFMSynthesis(fcVecMembrane, B, env, fs);

% EXAMPLE 3: feedback FM with pitch glide with center frequencies = membrane
% modal frequencies (YES)
% cool for bass drum, kettledrum
[yFBFMMemb3, yFBFMMembMat3] = feedbackFMSynthesis(fVecMembrane, BVec, env, fs);

% EXAMPLE 4: feedback FM with pitch glide with sounding frequencies = membrane
% modal frequencies
BVecLin = linspace(0.9, 0.905, N);
if strcmp(outDir, 'audioExamples/membrane/o-daiko/')==1 || ... 
            strcmp(outDir, 'audioExamples/membrane/snare/')==1 
    [yFBFMMemb4, yFBFMMembMat4] = feedbackFMSynthesis(fVecMembrane, BVecLin, env, fs);
else
    [yFBFMMemb4, yFBFMMembMat4] = feedbackFMSynthesis(fcVecMembrane, BVecLin, env, fs);
end

%% STRETCHED APF

b0 = (sqrt(1-B^2) - 1)/B;
fVecMembrane_w0 = fVecMembrane*sqrt(1 - B^2);

% EXAMPLE 1: wc = membrane modal frequencies
[ySAPFMemb1, ySAPFMembMat1] = stretchedAPFSynthesis(fVecMembrane_w0, b0, env, fs, [], 'none');

% EXAMPLE 2: w0 = membrane modal frequencies (YES)
% cool for bass drum
[ySAPFMemb2, ySAPFMembMat2] = stretchedAPFSynthesis(fVecMembrane, b0, env, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = membrane modal
% frequencies (YES)
if strcmp(outDir, 'audioExamples/membrane/o-daiko/')==1 || ...
            strcmp(outDir, 'audioExamples/membrane/snare/')==1 
    [ySAPFMemb3, ySAPFMembMat3] = stretchedAPFSynthesis(fVecMembrane_w0, b0, env, fs, fVecMembrane_w0 - 20*ones(1, Nf), 'linear');
else
    [ySAPFMemb3, ySAPFMembMat3] = stretchedAPFSynthesis(fVecMembrane, b0, env, fs, fVecMembrane - 20*ones(1, Nf), 'linear');
end
    
% EXAMPLE 4: feedback FM with pitch glide with w0 = membrane modal
% frequencies (YES)
fVecMembrane_w0End = fcVecMembrane .* sqrt(1 - BVecLin(end).^2);
fVecMembrane_w0EndK = fVecMembrane_w0 .* sqrt(1 - BVecLin(end).^2);
if strcmp(outDir, 'audioExamples/membrane/o-daiko/')==1 || ...
        strcmp(outDir, 'audioExamples/membrane/snare/')==1 
    [ySAPFMemb4, ySAPFMembMat4] = stretchedAPFSynthesis(fVecMembrane_w0, b0, env, fs, fVecMembrane_w0EndK, 'linear');
else
    [ySAPFMemb4, ySAPFMembMat4] = stretchedAPFSynthesis(fVecMembrane, b0, env, fs, fVecMembrane_w0End, 'linear');
end

%% TIME-VARYING APF

% EXAMPLE 1: time-varying APF center/sounding frequencies = membrane modal 
% frequencies, fixed APF parameters (YES!)
% cool for bass drum
[yTVAPFMemb1, yTVAPFMembMat1, TVAPFParams1] = TVAPFSynthesis(fVecMembrane, env, TVAPFParams, 0, fs);

% EXAMPLE 2: time-varying APF center/sounding frequencies = mesh modal 
% frequencies, randomized APF parameters - this sounds different everytime
% GOOD FOR SNARE
[yTVAPFMemb2, yTVAPFMembMat2, TVAPFParams2] = TVAPFSynthesis(fVecMembrane, env, TVAPFParams, 1, fs);

%% STRETCHED APF & TIME-VARYING APF vs. MODAL SYNTHESIS
% all time-varying APF calls here use randomized parameters

b0 = (sqrt(1-B^2) - 1)/B;
fVecMembrane_w0 = fVecMembrane*sqrt(1 - B^2);

% EXAMPLE 1: wc = mesh modal frequencies
[ySTVAPFMemb1, ySTVAPFMembMat1, ~] = stretchedAPFAndTVAPFSynthesis(fVecMembrane, b0, env, TVAPFParams, 0, fs, [], 'none');

% EXAMPLE 2: w0 = mesh modal frequencies
[ySTVAPFMemb2, ySTVAPFMembMat2, ~] = stretchedAPFAndTVAPFSynthesis(fVecMembrane_w0, b0, env, TVAPFParams, 0, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = mesh modal frequencies
[ySTVAPFMemb3, ySTVAPFMembMat3, ~] = stretchedAPFAndTVAPFSynthesis(fVecMembrane, b0, env, TVAPFParams, 0, fs, fVecMembrane - 20*ones(1, Nf), 'linear');

% EXAMPLE 4: feedback FM with pitch glide with w0 = mesh modal frequencies
% cool with the snare
[ySTVAPFMemb4, ySTVAPFMembMat4, ~] = stretchedAPFAndTVAPFSynthesis(fVecMembrane, b0, env, TVAPFParams, 0, fs, fVecMembrane_w0End, 'linear');

%% writing to file

if writeAudioFiles == 1
    
    if ~exist(outDir)
        mkdir(outDir)
    end
    
    audiowrite([outDir 'yMS.wav'], scaleForSavingAudio(real(yMS)), fs);

    audiowrite([outDir 'yFBFMMemb1.wav'], scaleForSavingAudio(real(yFBFMMemb1)), fs);
    audiowrite([outDir 'yFBFMMemb2.wav'], scaleForSavingAudio(real(yFBFMMemb2)), fs);
    audiowrite([outDir 'yFBFMMemb3.wav'], scaleForSavingAudio(real(yFBFMMemb3)), fs);
    audiowrite([outDir 'yFBFMMemb4.wav'], scaleForSavingAudio(real(yFBFMMemb4)), fs);
    
    audiowrite([outDir 'ySAPFMemb1.wav'], scaleForSavingAudio(real(ySAPFMemb1)), fs);
    audiowrite([outDir 'ySAPFMemb2.wav'], scaleForSavingAudio(real(ySAPFMemb2)), fs);
    audiowrite([outDir 'ySAPFMemb3.wav'], scaleForSavingAudio(real(ySAPFMemb3)), fs);
    audiowrite([outDir 'ySAPFMemb4.wav'], scaleForSavingAudio(real(ySAPFMemb4)), fs);
    
end