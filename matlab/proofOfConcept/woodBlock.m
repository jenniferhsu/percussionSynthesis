% This file creates sounds using the modal frequencies found in a wood
% block sample that I got from freesound.org.
%% input parameters 

% general
fs = 44100;
dur = 1;
plotSpectrograms = 0;
writeAudioFiles = 1;

% feedback FM
B = 0.9;    % feedback coefficient
g = 0.9999; % pitch glide coefficient
gWood = 0.9993;

% time-varying APF
TVAPFParams.M = fs/40;
TVAPFParams.f_m = 100;
TVAPFParams.f_b = fs/16;

%% derived parameters
N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient

%% analyze a wood block recording to get the modal frequencies
[wb, fs2] = audioread('audioExamples/woodBlocks/218461__thomasjaunism__wood-block.wav');
wb = wb(:,1);

Nwb = length(wb);
Nfft = 2^nextpow2(Nwb);
faxis = (fs2/2) * linspace(0, 1, Nfft/2+1);

WB = fft(wb, Nfft);
WBPos = WB(1:Nfft/2+1);
WBPosdB = 20*log10(abs(WBPos)/max(abs(WBPos)));


[pks, locs] = findpeaks(WBPosdB, 'minpeakheight', -60, 'minpeakdist', 80);
[spks, sinds] = sort(pks, 'descend');
spks = spks(1:11);
slocs = locs(sinds(1:11));

plot(faxis, WBPosdB);
hold on
for i=1:length(spks)
    plot(faxis(slocs(i)), spks(i), 'r*');
end

fVecWB = faxis(slocs);
fVecWB = sort(fVecWB, 'ascend');
fVecWB = fVecWB(2:end); % drop the lowest one at like 5Hz

Nf = length(fVecWB);


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
eWood = gWood.^(linspace(0, N, N)); % exponential decay for woodblock

envWood = zeros(Nf, N);
for i=1:Nf
    envWood(i,:) = aStart(i) * eWood;
end

%% modal synthesis with exponential decay using the modes
yMS = zeros(1, N);
nVec = 0:T:(dur-T);

for i=1:Nf
    f = fVecWB(i);
    yMS = yMS + (exp(1j*2*pi*f*nVec) .* envWood(i,:));
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yMS), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
end

%% FEEDBACK FM

% EXAMPLE 1: feedback FM carrier frequencies = membrane modal frequencies
% (YES)
[yFBFMWB1, yFBFMWBMat1] = feedbackFMSynthesis(fVecWB, B, envWood, fs);

% EXAMPLE 3: feedback FM with pitch glide with carrier frequencies = membrane
% modal frequencies (YES)
% cool for bass drum
[yFBFMWB2, yFBFMWBMat2] = feedbackFMSynthesis(fVecWB, BVec, envWood, fs);



%% STRETCHED APF

B = 0.9;
b0 = (sqrt(1-B^2) - 1)/B;
fVecWB_w0 = fVecWB*sqrt(1 - B^2);

% EXAMPLE 1: w0 = membrane modal frequencies
[ySAPFWB1, ySAPFWBMat1] = stretchedAPFSynthesis(fVecWB, b0, envWood, fs, [], 'none');

% EXAMPLE 2: feedback FM with FBFM pitch glide with wc = membrane modal
% frequencies 
%[ySAPFWB2, ySAPFWBMat2] = stretchedAPFSynthesis(fVecWB_w0, b0, envWood, fs, 2*fVecWB_w0, 'linear');
%[ySAPFWB2, ySAPFWBMat2] = stretchedAPFSynthesis(fVecWB, b0, envWood, fs, 2*fVecWB, 'linear');
[ySAPFWB2, ySAPFWBMat2] = stretchedAPFSynthesis(fVecWB, b0, envWood, fs, 4*fVecWB, 'exp');

% EXAMPLE 3: feedback FM with pitch glide with w0 = membrane modal
% frequencies (YES)
[ySAPFWB3, ySAPFWBMat3] = stretchedAPFSynthesis(fVecWB, b0, envWood, fs, fVecWB*1.2, 'linear');


%% TIME-VARYING APF

% EXAMPLE 1: time-varying APF center/sounding frequencies = membrane modal 
% frequencies, fixed APF parameters (YES!)
% cool for bass drum
[yTVAPFWB1, yTVAPFWBMat1, TVAPFParams1] = TVAPFSynthesis(fVecWB, envWood, TVAPFParams, 0, fs);

% EXAMPLE 2: time-varying APF center/sounding frequencies = mesh modal 
% frequencies, randomized APF parameters - this sounds different everytime
% GOOD FOR SNARE
[yTVAPFWB2, yTVAPFWBMat2, TVAPFParams2] = TVAPFSynthesis(fVecWB, envWood, TVAPFParams, 1, fs);

%% STRETCHED APF & TIME-VARYING APF vs. MODAL SYNTHESIS
% all time-varying APF calls here use randomized parameters

b0 = (sqrt(1-B^2) - 1)/B;
fVecWB_wc = fVecWB*sqrt(1 - B^2);

% EXAMPLE 1: wc = mesh modal frequencies
[ySTVAPFWB1, ySTVAPFWBMat1, ~] = stretchedAPFAndTVAPFSynthesis(fVecWB, b0, envWood, TVAPFParams, 0, fs, [], 'none');

% EXAMPLE 2: w0 = mesh modal frequencies
[ySTVAPFWB2, ySTVAPFWBMat2, ~] = stretchedAPFAndTVAPFSynthesis(fVecWB_wc, b0, envWood, TVAPFParams, 0, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = mesh modal frequencies
[ySTVAPFWB3, ySTVAPFWBMat3, ~] = stretchedAPFAndTVAPFSynthesis(fVecWB, b0, envWood, TVAPFParams, 0, fs, fVecWB, 'fbfm');

% EXAMPLE 4: feedback FM with pitch glide with w0 = mesh modal frequencies
% cool with the snare
[ySTVAPFWB4, ySTVAPFWBMat5, ~] = stretchedAPFAndTVAPFSynthesis(fVecWB, b0, envWood, TVAPFParams, 0, fs, fVecWB*1.2, 'linear');

%% Write to file


if writeAudioFiles == 1
    
   outputDir = 'audioExamples/woodBlocks/';
    if ~exist(outputDir)
        mkdir(outputDir);
    end
    
    % steel plate modal synthesis
    audiowrite([outputDir 'yMS.wav'], scaleForSavingAudio(real(yMS)), fs)
    
    % feedback FM vs. steel plate modal synthesis
    audiowrite([outputDir 'yFBFMWB1.wav'], scaleForSavingAudio(real(yFBFMWB1)), fs)
    audiowrite([outputDir 'yFBFMWB2.wav'], scaleForSavingAudio(real(yFBFMWB2)), fs)
    
    % stretched APF vs. steel plate modal synthesis
    audiowrite([outputDir 'ySAPFWB1.wav'], scaleForSavingAudio(real(ySAPFWB1)), fs)
    audiowrite([outputDir 'ySAPFWB2.wav'], scaleForSavingAudio(real(ySAPFWB2)), fs)
    audiowrite([outputDir 'ySAPFWB3.wav'], scaleForSavingAudio(real(ySAPFWB3)), fs)
    
    % time-varying APF vs. steel plate modal synthesis
    audiowrite([outputDir 'yTVAPFWB1.wav'], scaleForSavingAudio(real(yTVAPFWB1)), fs)
    audiowrite([outputDir 'yTVAPFWB2.wav'], scaleForSavingAudio(real(yTVAPFWB2)), fs)
    
    % stretched APF and time-varying APF vs. steel plate modal synthesis
    audiowrite([outputDir 'ySTVAPFWB1.wav'], scaleForSavingAudio(real(ySTVAPFWB1)), fs)
    audiowrite([outputDir 'ySTVAPFWB2.wav'], scaleForSavingAudio(real(ySTVAPFWB2)), fs)
    audiowrite([outputDir 'ySTVAPFWB3.wav'], scaleForSavingAudio(real(ySTVAPFWB3)), fs)
    audiowrite([outputDir 'ySTVAPFWB4.wav'], scaleForSavingAudio(real(ySTVAPFWB4)), fs)
    
end
