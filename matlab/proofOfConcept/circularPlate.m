% This file creates sounds using circular plate modal frequencies

%% input parameters 

% general
fs = 44100;
dur = 1;
plotSpectrograms = 0;
writeAudioFiles = 1;

% feedback FM
B = 0.99;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

% time-varying APF
TVAPFParams.M = fs/40;
TVAPFParams.f_m = 100;
TVAPFParams.f_b = fs/16;

%% derived parameters
N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient

%% circular plate with a clamped edge
% from Science of Percussion Instruments (page 80)

h = 0.005;
a = 0.09;

E = 2*10^11;
v = 0.3;
rho = 7860;
cL = sqrt(E/(rho*(1 - v^2)));

%f01 = 0.4694*cL*h/a^2;
%f01 = 2000;
%fVecCP = f01 * [1, 2.08, 3.41, 5, 6.82, 3.89, 5.95, 8.28, 10.87, 13.71, 8.72, 11.75, 15.06, 18.63, 22.47];
% SIMPLY-SUPPORTED CIRCULAR PLATE

f01 = 0.2287 * cL * h / a^2;
fVecCP = f01 * [1, 2.80, 5.15, 5.98, 9.75, 14.09, 14.91, 20.66, 26.99];

Nf = length(fVecCP);

% center frequencies for feedback FM if sounding frequencies = membrane
% modal freqeuncies
fcVecCP = fVecCP./(sqrt(1-B^2));

%% set up decay envelopes
% want the lowest modal frequencies to have a larger starting amplitude 
% than the highest modal frequencies
aStart = zeros(1, Nf);
aStart(1) = 1;              % highest amplitude for lowest modal freq
aStart(Nf) = 0.5;          % lowest amplitude for highest modal freq
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
    f = fVecCP(i);
    yMS = yMS + (exp(1j*2*pi*f*nVec) .* env(i,:));
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yMS), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
end


%% FEEDBACK FM

% EXAMPLE 1: feedback FM center frequencies = membrane modal frequencies
% (YES)
[yFBFMCP1, yFBFMCPMat1] = feedbackFMSynthesis(fVecCP, B, env, fs);

% EXAMPLE 2: feedback FM sounding frequencies = membrane modal frequencies
[yFBFMCP2, yFBFMCPMat2] = feedbackFMSynthesis(fcVecCP, B, env, fs);

% EXAMPLE 3: feedback FM with pitch glide with center frequencies = membrane
% modal frequencies (YES)
% cool for bass drum
BVec = linspace(0.99, 0.989, N);
[yFBFMCP3, yFBFMCPMat3] = feedbackFMSynthesis(fcVecCP, BVec, env, fs);

% EXAMPLE 4: feedback FM with pitch glide with sounding frequencies = membrane
% modal frequencies
BVec = linspace(0.91, 0.9, N);
[yFBFMCP4, yFBFMCPMat4] = feedbackFMSynthesis(fcVecCP, BVec, env, fs);


%% STRETCHED APF

b0 = (sqrt(1-B^2) - 1)/B;
fVecCP_wc = fVecCP*sqrt(1 - B^2);

% EXAMPLE 1: wc = membrane modal frequencies
[ySAPFCP1, ySAPFCPMat1] = stretchedAPFSynthesis(fVecCP, b0, env, fs, [], 'none');

% EXAMPLE 2: w0 = membrane modal frequencies (YES)
% cool for bass drum
[ySAPFCP2, ySAPFCPMat2] = stretchedAPFSynthesis(fcVecCP, b0, env, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = membrane modal
% frequencies 
BVec = linspace(0.99, 0.989, N);
fcVecCPStart = fcVecCP .* sqrt(1 - BVec(1)^2);
fcVecCPEnd = fcVecCP .* sqrt(1 - BVec(end)^2);
[ySAPFCP3, ySAPFCPMat3] = stretchedAPFSynthesis(fcVecCPStart, b0, env, fs, fcVecCPEnd, 'linear');

% EXAMPLE 4: feedback FM with pitch glide with w0 = membrane modal
% frequencies (YES)
BVec = linspace(0.989, 0.99, N);
fcVecCPStart = fcVecCP .* sqrt(1 - BVec(1)^2);
fcVecCPEnd = fcVecCP .* sqrt(1 - BVec(end)^2);
[ySAPFCP4, ySAPFCPMat4] = stretchedAPFSynthesis(fcVecCPStart, b0, env, fs, fcVecCPEnd, 'linear');

%% TIME-VARYING APF

% EXAMPLE 1: time-varying APF center/sounding frequencies = membrane modal 
% frequencies, fixed APF parameters (YES!)
% cool for bass drum
[yTVAPFCP1, yTVAPFCPMat1, TVAPFParams1] = TVAPFSynthesis(fVecCP, env, TVAPFParams, 0, fs);

% EXAMPLE 2: time-varying APF center/sounding frequencies = mesh modal 
% frequencies, randomized APF parameters - this sounds different everytime
% GOOD FOR SNARE
[yTVAPFCP2, yTVAPFCPMat2, TVAPFParams2] = TVAPFSynthesis(fVecCP, env, TVAPFParams, 1, fs);

%% STRETCHED APF & TIME-VARYING APF vs. MODAL SYNTHESIS
% all time-varying APF calls here use randomized parameters

b0 = (sqrt(1-B^2) - 1)/B;
fVecCP_wc = fVecCP*sqrt(1 - B^2);

% EXAMPLE 1: wc = mesh modal frequencies
[ySTVAPFCP1, ySTVAPFCPMat1, ~] = stretchedAPFAndTVAPFSynthesis(fVecCP, b0, env, TVAPFParams, 0, fs, [], 'none');

% EXAMPLE 2: w0 = mesh modal frequencies
[ySTVAPFCP2, ySTVAPFCPMat2, ~] = stretchedAPFAndTVAPFSynthesis(fVecCP_wc, b0, env, TVAPFParams, 0, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = mesh modal frequencies
[ySTVAPFCP3, ySTVAPFCPMat3, ~] = stretchedAPFAndTVAPFSynthesis(fVecCP, b0, env, TVAPFParams, 0, fs, fVecCP, 'fbfm');

% EXAMPLE 4: feedback FM with pitch glide with w0 = mesh modal frequencies
% cool with the snare
[ySTVAPFCP4, ySTVAPFCPMat4, ~] = stretchedAPFAndTVAPFSynthesis(fVecCP, b0, env, TVAPFParams, 0, fs, fVecCP*1.2, 'linear');

%% writing to file

if writeAudioFiles == 1
    
   outputDir = 'audioExamples/circularPlate/';
    if ~exist(outputDir)
        mkdir(outputDir);
    end
    
    % steel plate modal synthesis
    audiowrite([outputDir 'yMS.wav'], scaleForSavingAudio(real(yMS)), fs)
    
    % feedback FM vs. steel plate modal synthesis
    audiowrite([outputDir 'yFBFMCP1.wav'], scaleForSavingAudio(real(yFBFMCP1)), fs)
    audiowrite([outputDir 'yFBFMCP2.wav'], scaleForSavingAudio(real(yFBFMCP2)), fs)
    audiowrite([outputDir 'yFBFMCP3.wav'], scaleForSavingAudio(real(yFBFMCP3)), fs)
    audiowrite([outputDir 'yFBFMCP4.wav'], scaleForSavingAudio(real(yFBFMCP4)), fs)
    
    % stretched APF vs. steel plate modal synthesis
    audiowrite([outputDir 'ySAPFCP1.wav'], scaleForSavingAudio(real(ySAPFCP1)), fs)
    audiowrite([outputDir 'ySAPFCP2.wav'], scaleForSavingAudio(real(ySAPFCP2)), fs)
    audiowrite([outputDir 'ySAPFCP3.wav'], scaleForSavingAudio(real(ySAPFCP3)), fs)
    audiowrite([outputDir 'ySAPFCP4.wav'], scaleForSavingAudio(real(ySAPFCP4)), fs)
    
    % time-varying APF vs. steel plate modal synthesis
    audiowrite([outputDir 'yTVAPFCP1.wav'], scaleForSavingAudio(real(yTVAPFCP1)), fs)
    audiowrite([outputDir 'yTVAPFCP2.wav'], scaleForSavingAudio(real(yTVAPFCP2)), fs)
    
    % stretched APF and time-varying APF vs. steel plate modal synthesis
    audiowrite([outputDir 'ySTVAPFCP1.wav'], scaleForSavingAudio(real(ySTVAPFCP1)), fs)
    audiowrite([outputDir 'ySTVAPFCP2.wav'], scaleForSavingAudio(real(ySTVAPFCP2)), fs)
    audiowrite([outputDir 'ySTVAPFCP3.wav'], scaleForSavingAudio(real(ySTVAPFCP3)), fs)
    audiowrite([outputDir 'ySTVAPFCP4.wav'], scaleForSavingAudio(real(ySTVAPFCP4)), fs)
    
end
