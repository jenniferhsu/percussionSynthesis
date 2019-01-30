% This file creates sounds using the rectangular plate equation
% we make a wood block sound, make it sound more like glass, and then like
% a metallic sheet
%% input parameters 

% general
fs = 44100;
dur = 1;
plotSpectrograms = 0;
writeAudioFiles = 1;

% feedback FM
B = 0.9;    % feedback coefficient
g = 0.9999; % pitch glide coefficient
gWood = 0.997;

% time-varying APF
TVAPFParams.M = fs/40;
TVAPFParams.f_m = 100;
TVAPFParams.f_b = fs/16;

%% derived parameters
N = fs*dur;
T = 1/fs;
envWood = gWood.^(linspace(0, N, N)); % exponential decay for woodblock
env = g.^(linspace(0, N, N));   % exponential decay
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient

%% rectangular plate with a simply supported edge
% from Science of Percussion Instruments (page 81)

f01 = 180;
Lx = .9;
Ly = .9;
%Ly = 1.6;
mm = 5;
nn = 5;
fVec = zeros(1, (mm+1)*(nn+1));
i = 1;
for m=0:mm
    for n=0:nn
        fVec(i) = ((m+1)/Lx)^2 + ((n+1)/Ly)^2; 
        i = i+1;
    end
end
fVec = fVec/fVec(1);
fVec = unique(fVec);

fVecRP = f01 * fVec;

Nf = length(fVecRP);

% center frequencies for feedback FM if sounding frequencies = membrane
% modal freqeuncies
fcVecRP = fVecRP./(sqrt(1-B^2));

%% modal synthesis with exponential decay using the modes
yMS = zeros(1, N);
nVec = 0:T:(dur-T);

for i=1:Nf
    f = fVecRP(i);
    yMS = yMS + (exp(1j*2*pi*f*nVec) .* envWood);
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yMS), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
end

%% FEEDBACK FM

% EXAMPLE 1: feedback FM center frequencies = membrane modal frequencies
% (YES)
[yFBFMRP1, yFBFMRPMat1] = feedbackFMSynthesis(fVecRP, B, envWood, fs);

% EXAMPLE 2: feedback FM sounding frequencies = membrane modal frequencies
[yFBFMRP2, yFBFMRPMat2] = feedbackFMSynthesis(fcVecRP, B, envWood, fs);

% EXAMPLE 3: feedback FM with pitch glide with center frequencies = membrane
% modal frequencies (YES)
% cool for bass drum
[yFBFMRP3, yFBFMRPMat3] = feedbackFMSynthesis(fVecRP, BVec, envWood, fs);

% EXAMPLE 4: feedback FM with pitch glide with sounding frequencies = membrane
% modal frequencies
[yFBFMRP4, yFBFMRPMat4] = feedbackFMSynthesis(fcVecRP, BVec, envWood, fs);


%% STRETCHED APF

b0 = (sqrt(1-B^2) - 1)/B;
fVecRP_wc = fVecRP*sqrt(1 - B^2);

% EXAMPLE 1: wc = membrane modal frequencies
[ySAPFRP1, ySAPFRPMat1] = stretchedAPFSynthesis(fVecRP, b0, envWood, fs, [], 'none');

% EXAMPLE 2: w0 = membrane modal frequencies (YES)
% cool for bass drum
[ySAPFRP2, ySAPFRPMat2] = stretchedAPFSynthesis(fVecRP_wc, b0, envWood, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = membrane modal
% frequencies (YES)
[ySAPFRP3, ySAPFRPMat3] = stretchedAPFSynthesis(fVecRP, b0, envWood, fs, fVecRP*2, 'fbfm');

% EXAMPLE 4: feedback FM with pitch glide with w0 = membrane modal
% frequencies (YES)
[ySAPFRP4, ySAPFRPMat4] = stretchedAPFSynthesis(fVecRP, b0, envWood, fs, fVecRP*1.2, 'linear');


%% TIME-VARYING APF

% EXAMPLE 1: time-varying APF center/sounding frequencies = membrane modal 
% frequencies, fixed APF parameters (YES!)
% cool for bass drum
[yTVAPFRP1, yTVAPFRPMat1, TVAPFParams1] = TVAPFSynthesis(fVecRP, envWood, TVAPFParams, 0, fs);

% EXAMPLE 2: time-varying APF center/sounding frequencies = mesh modal 
% frequencies, randomized APF parameters - this sounds different everytime
% GOOD FOR SNARE
[yTVAPFRP2, yTVAPFRPMat2, TVAPFParams2] = TVAPFSynthesis(fVecRP, envWood, TVAPFParams, 1, fs);

%% STRETCHED APF & TIME-VARYING APF vs. MODAL SYNTHESIS
% all time-varying APF calls here use randomized parameters

b0 = (sqrt(1-B^2) - 1)/B;
fVecRP_wc = fVecRP*sqrt(1 - B^2);

% EXAMPLE 1: wc = mesh modal frequencies
[ySTVAPFRP1, ySTVAPFRPMat1, ~] = stretchedAPFAndTVAPFSynthesis(fVecRP, b0, envWood, TVAPFParams, 0, fs, [], 'none');

% EXAMPLE 2: w0 = mesh modal frequencies
[ySTVAPFRP2, ySTVAPFRPMat2, ~] = stretchedAPFAndTVAPFSynthesis(fVecRP_wc, b0, envWood, TVAPFParams, 0, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = mesh modal frequencies
[ySTVAPFRP3, ySTVAPFRPMat3, ~] = stretchedAPFAndTVAPFSynthesis(fVecRP, b0, envWood, TVAPFParams, 0, fs, fVecRP, 'fbfm');

% EXAMPLE 4: feedback FM with pitch glide with w0 = mesh modal frequencies
% cool with the snare
[ySTVAPFRP4, ySTVAPFRPMat5, ~] = stretchedAPFAndTVAPFSynthesis(fVecRP, b0, envWood, TVAPFParams, 0, fs, fVecRP*1.2, 'linear');

%% Write to file


if writeAudioFiles == 1
    
   outputDir = 'audioExamples/woodBlocks/';
    if ~exist(outputDir)
        mkdir(outputDir);
    end
    
    % steel plate modal synthesis
    audiowrite([outputDir 'yMS.wav'], scaleForSavingAudio(real(yMS)), fs)
    
    % feedback FM vs. steel plate modal synthesis
    audiowrite([outputDir 'yFBFMRP1.wav'], scaleForSavingAudio(real(yFBFMRP1)), fs)
    audiowrite([outputDir 'yFBFMRP2.wav'], scaleForSavingAudio(real(yFBFMRP2)), fs)
    audiowrite([outputDir 'yFBFMRP3.wav'], scaleForSavingAudio(real(yFBFMRP3)), fs)
    audiowrite([outputDir 'yFBFMRP4.wav'], scaleForSavingAudio(real(yFBFMRP4)), fs)
    
    % stretched APF vs. steel plate modal synthesis
    audiowrite([outputDir 'ySAPFRP1.wav'], scaleForSavingAudio(real(ySAPFRP1)), fs)
    audiowrite([outputDir 'ySAPFRP2.wav'], scaleForSavingAudio(real(ySAPFRP2)), fs)
    audiowrite([outputDir 'ySAPFRP3.wav'], scaleForSavingAudio(real(ySAPFRP3)), fs)
    audiowrite([outputDir 'ySAPFRP4.wav'], scaleForSavingAudio(real(ySAPFRP4)), fs)
    
    % time-varying APF vs. steel plate modal synthesis
    audiowrite([outputDir 'yTVAPFRP1.wav'], scaleForSavingAudio(real(yTVAPFRP1)), fs)
    audiowrite([outputDir 'yTVAPFRP2.wav'], scaleForSavingAudio(real(yTVAPFRP2)), fs)
    
    % stretched APF and time-varying APF vs. steel plate modal synthesis
    audiowrite([outputDir 'ySTVAPFRP1.wav'], scaleForSavingAudio(real(ySTVAPFRP1)), fs)
    audiowrite([outputDir 'ySTVAPFRP2.wav'], scaleForSavingAudio(real(ySTVAPFRP2)), fs)
    audiowrite([outputDir 'ySTVAPFRP3.wav'], scaleForSavingAudio(real(ySTVAPFRP3)), fs)
    audiowrite([outputDir 'ySTVAPFRP4.wav'], scaleForSavingAudio(real(ySTVAPFRP4)), fs)
    
end
