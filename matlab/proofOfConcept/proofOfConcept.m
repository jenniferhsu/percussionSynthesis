% proofOfConcept.m
% 
% This script is a proof of concept that the feedback FM and time-varying
% allpass filters technique can be used to generate interesting percussive
% sounds. Sounds are generated for
%   a) the 3x3 DWG mesh with a reflection coefficient at the boundaries
%   b) the 3x3 DWG mesh with a one-zero lowpass filter at the boundaries
%   c) feedback FM with center frequency wc equal to the modal frequencies
%       from (a)
%   d) feedback FM with center frequency wc equal to the modal frequencies
%       from (a)
%   e) feedback FM with a pitch glide with center frequency wc equal to the 
%       modal frequencies from (a)
%   f) feedback FM with a pitch glide sounding frequency wc equal to the 
%       modal frequencies from (a)
%   g) modal synthesis using the frequencies found from the 3x3 mode steel 
%       plate von Karman equations
%   h) feedback FM with center frequency wc equal to the modal frequencies
%       from (g)
%   i) feedback FM with sounding frequency wc equal to the modal frequencies
%       from (g)
%   j) feedback FM with a pitch glide with center frequency wc equal to the 
%       modal frequencies from (g)
%   k) feedback FM with a pitch glide sounding frequency wc equal to the 
%       modal frequencies from (g)
%
% author: Jennifer Hsu
% date: 1/11/2019

addpath(genpath('../basicMeshFunctions'));
addpath(genpath('../modalFrequencyEquations'));

%% input parameters 

% general
fs = 44100;
dur = 1;
plotSpectrograms = 0;
writeAudioFiles = 0;

% mesh
ex = 1;
Njr = 3;
Njc = 3;
inX = 3;
inY = 3;
outX = 3;
outY = 3;
refl = -0.999;
BLPF = -[0.5 0.5];
ALPF = 1;
plotOn = 0;

% feedback FM
f1 = 5500;  % sounding frequency of a 3x3 digital waveguide mesh, Hz
B = 0.9;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

% time-varying APF
TVAPFParams.M = fs/40;
TVAPFParams.f_m = 100;
TVAPFParams.f_b = fs/16;


%% derived parameters
N = fs*dur;
T = 1/fs;
%env = g.^(linspace(0, N, N));   % exponential decay
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient


%% MESH

% 3x3 mesh using the reflection coefficient
yMeshReflCoeff = mesh2DRectReflCoeff(ex, Njr, Njc, N, fs, inX, inY, outX, outY, refl, plotOn);

% 3x3 mesh using a one-zero low pass filter
yMeshFilter = mesh2DRectFilterFast(ex, Njr, Njc, N, fs, inX, inY, outX, outY, BLPF, ALPF, plotOn);

if plotSpectrograms == 1
    figure
    spectrogram(real(yMeshReflCoeff), hann(256), 128, 1024, fs, 'yaxis');
    title('3x3 mesh with reflection coefficient spectrogram')

    figure
    spectrogram(real(yMeshFilter), hann(256), 128, 1024, fs, 'yaxis');
    title('3x3 mesh with one-zero LPF spectrogram')
end

% Fourier analysis of 3x3 mesh with a reflection coefficient to pick out
% the modal frequencies. The one with filters at the boundaries doesn't
% last long enough for this kind of analysis.
Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);
Y = fft(yMeshReflCoeff, Nfft);
YPos = Y(1:Nfft/2+1);
[pks, locs] = findpeaks(abs(YPos), 'minpeakheight', 50, 'minpeakdist', 5);

figure
plot(faxis, abs(YPos))
hold on
for l=1:length(locs)
    plot(faxis(locs(l)), pks(l), 'r*');
end

% mesh modal frequencies
fVecMesh = faxis(locs);
Nf = length(fVecMesh);

% center frequencies for feedback FM if sounding frequencies = mesh modal
% frequencies
fcVecMesh = fVecMesh./(sqrt(1-B^2));

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


%% FEEDBACK FM vs. MESH

% EXAMPLE 1: feedback FM center frequencies = mesh modal frequencies
[yFBFMMesh1, yFBFMMeshMat1] = feedbackFMSynthesis(fVecMesh, B, env, fs);

% EXAMPLE 2: feedback FM sounding frequencies = mesh modal frequencies
[yFBFMMesh2, yFBFMMeshMat2] = feedbackFMSynthesis(fcVecMesh, B, env, fs);

% EXAMPLE 3: feedback FM with pitch glide with center frequencies = mesh
% modal frequencies
[yFBFMMesh3, yFBFMMeshMat3] = feedbackFMSynthesis(fVecMesh, BVec, env, fs);

% EXAMPLE 4: feedback FM with pitch glide with sounding frequencies = mesh
% modal frequencies
[yFBFMMesh4, yFBFMMeshMat4] = feedbackFMSynthesis(fcVecMesh, BVec, env, fs);

if plotSpectrograms == 1
    figure
    spectrogram(real(yFBFMMesh1), hann(256), 128, 1024, fs, 'yaxis');
    title('feedback FM synthesis using 3x3 mesh modal frequencies spectrogram 1')

    figure
    spectrogram(real(yFBFMMesh2), hann(256), 128, 1024, fs, 'yaxis');
    title('feedback FM synthesis using 3x3 mesh modal frequencies spectrogram 2')

    figure
    spectrogram(real(yFBFMMesh3), hann(256), 128, 1024, fs, 'yaxis');
    title('feedback FM synthesis using 3x3 mesh modal frequencies spectrogram 3')

    figure
    spectrogram(real(yFBFMMesh4), hann(256), 128, 1024, fs, 'yaxis');
    title('feedback FM synthesis using 3x3 mesh modal frequencies spectrogram 4') 
end

%% STRETCHED APF vs. MESH

b0 = (sqrt(1-B^2) - 1)/B;
fVecMesh_wc = fVecMesh*sqrt(1 - B^2);

% EXAMPLE 1: wc = mesh modal frequencies
[ySAPFMesh1, ySAPFMeshMat1] = stretchedAPFSynthesis(fVecMesh, b0, env, fs, [], 'none');

% EXAMPLE 2: w0 = mesh modal frequencies
[ySAPFMesh2, ySAPFMeshMat2] = stretchedAPFSynthesis(fVecMesh_wc, b0, env, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = mesh modal frequencies
[ySAPFMesh3, ySAPFMeshMat3] = stretchedAPFSynthesis(fVecMesh, b0, env, fs, fVecMesh, 'fbfm');

% EXAMPLE 4: feedback FM with pitch glide with w0 = mesh modal frequencies
[ySAPFMesh4, ySAPFMeshMat4] = stretchedAPFSynthesis(fVecMesh_wc, b0, env, fs, fVecMesh_wc, 'fbfm');

if plotSpectrograms == 1
    figure
    spectrogram(real(ySAPFMesh1), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using 3x3 mesh modal frequencies spectrogram 1')

    figure
    spectrogram(real(ySAPFMesh2), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using 3x3 mesh modal frequencies spectrogram 2')

    figure
    spectrogram(real(ySAPFMesh3), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using 3x3 mesh modal frequencies spectrogram 3')

    figure
    spectrogram(real(ySAPFMesh4), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using 3x3 mesh modal frequencies spectrogram 4') 
end

%% TIME-VARYING APF vs. MESH

% EXAMPLE 1: time-varying APF center/sounding frequencies = mesh modal 
% frequencies, fixed APF parameters
[yTVAPFMesh1, yTVAPFMeshMat1, TVAPFParams1] = TVAPFSynthesis(fVecMesh, env, TVAPFParams, 0, fs);

% EXAMPLE 2: time-varying APF center/sounding frequencies = mesh modal 
% frequencies, randomized APF parameters - this sounds different everytime
[yTVAPFMesh2, yTVAPFMeshMat2, TVAPFParams2] = TVAPFSynthesis(fVecMesh, env, TVAPFParams, 1, fs);

if plotSpectrograms == 1
    figure
    spectrogram(real(yTVAPFMesh1), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying APF synthesis using 3x3 mesh modal frequencies spectrogram - fixed params')

    figure
    spectrogram(real(yTVAPFMesh2), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying APF synthesis using 3x3 mesh modal frequencies spectrogram - random params')
end

%% STRETCHED APF & TIME-VARYING APF vs. MESH
% all time-varying APF calls here use randomized parameters

b0 = (sqrt(1-B^2) - 1)/B;
fVecMesh_wc = fVecMesh*sqrt(1 - B^2);

% EXAMPLE 1: wc = mesh modal frequencies
[ySTVAPFMesh1, ySTVAPFMeshMat1, ~] = stretchedAPFAndTVAPFSynthesis(fVecMesh, b0, env, TVAPFParams, 1, fs, [], 'none');

% EXAMPLE 2: w0 = mesh modal frequencies
[ySTVAPFMesh2, ySTVAPFMeshMat2, ~] = stretchedAPFAndTVAPFSynthesis(fVecMesh_wc, b0, env, TVAPFParams, 1, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = mesh modal frequencies
[ySTVAPFMesh3, ySTVAPFMeshMat3, ~] = stretchedAPFAndTVAPFSynthesis(fVecMesh, b0, env, TVAPFParams, 1, fs, fVecMesh, 'fbfm');

% EXAMPLE 4: feedback FM with pitch glide with w0 = mesh modal frequencies
[ySTVAPFMesh4, ySTVAPFMeshMat4, ~] = stretchedAPFAndTVAPFSynthesis(fVecMesh, b0, env, TVAPFParams, 1, fs, fVecMesh_wc, 'fbfm');

if plotSpectrograms == 1
    figure
    spectrogram(real(ySTVAPFMesh1), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using 3x3 mesh modal frequencies spectrogram 1')

    figure
    spectrogram(real(ySTVAPFMesh2), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using 3x3 mesh modal frequencies spectrogram 2')

    figure
    spectrogram(real(ySTVAPFMesh3), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using 3x3 mesh modal frequencies spectrogram 3')

    figure
    spectrogram(real(ySTVAPFMesh4), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using 3x3 mesh modal frequencies spectrogram 4') 
end



%% MODAL SYNTHESIS

% generate modes for a 3x3 mode steel plate using the von Karman equations
% (fVec is used for feedback FM synthesis as well!)
modes = plateRectModes(Njr, Njc, 1, 1);
modes = modes(:);
modes = unique(modes);
fVecModal = modes*f1;
Nf = length(fVecModal);

% center frequencies for feedback FM if sounding frequencies = von Karman 
%modal frequencies
fcVecModal = fVecModal./(sqrt(1-B^2));

% basic modal/additive synthesis with exponential decay using the modes
yMS = zeros(1, N);
nVec = 0:T:(dur-T);

for i=1:Nf
    f = fVecModal(i);
    yMS = yMS + (exp(1j*2*pi*f*nVec) .* env);
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yMS), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
end

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


%% FEEDBACK FM vs. MODAL SYNTHESIS

% EXAMPLE 1: feedback FM center frequencies = mesh modal frequencies
[yFBFMModal1, yFBFMModalMat1] = feedbackFMSynthesis(fVecModal, B, env, fs);

% EXAMPLE 2: feedback FM sounding frequencies = mesh modal frequencies
[yFBFMModal2, yFBFMModalMat2] = feedbackFMSynthesis(fcVecModal, B, env, fs);

% EXAMPLE 3: feedback FM with pitch glide with center frequencies = mesh
% modal frequencies
[yFBFMModal3, yFBFMModalMat3] = feedbackFMSynthesis(fVecModal, BVec, env, fs);

% EXAMPLE 4: feedback FM with pitch glide with sounding frequencies = mesh
% modal frequencies
[yFBFMModal4, yFBFMModalMat4] = feedbackFMSynthesis(fcVecModal, BVec, env, fs);

if plotSpectrograms == 1
    figure
    spectrogram(real(yFBFMModal1), hann(256), 128, 1024, fs, 'yaxis');
    title('feedback FM synthesis using theoretical steel plate modal frequencies spectrogram 1')

    figure
    spectrogram(real(yFBFMModal2), hann(256), 128, 1024, fs, 'yaxis');
    title('feedback FM synthesis using theoretical steel plate modal frequencies spectrogram 2')

    figure
    spectrogram(real(yFBFMModal3), hann(256), 128, 1024, fs, 'yaxis');
    title('feedback FM synthesis using theoretical steel plate modal frequencies spectrogram 3')

    figure
    spectrogram(real(yFBFMModal4), hann(256), 128, 1024, fs, 'yaxis');
    title('feedback FM synthesis using theoretical steel plate modal frequencies spectrogram 4')
end

%% STRETCHED APF vs. MODAL SYNTHESIS

b0 = (sqrt(1-B^2) - 1)/B;
fVecModal_wc = fVecModal*sqrt(1 - B^2);

% EXAMPLE 1: wc = mesh modal frequencies
[ySAPFModal1, ySAPFModalMat1] = stretchedAPFSynthesis(fVecModal, b0, env, fs, [], 'none');

% EXAMPLE 2: w0 = mesh modal frequencies
[ySAPFModal2, ySAPFModalMat2] = stretchedAPFSynthesis(fVecModal_wc, b0, env, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = mesh modal frequencies
[ySAPFModal3, ySAPFModalMat3] = stretchedAPFSynthesis(fVecModal, b0, env, fs, fVecModal, 'fbfm');

% EXAMPLE 4: feedback FM with pitch glide with w0 = mesh modal frequencies
[ySAPFModal4, ySAPFModalMat4] = stretchedAPFSynthesis(fVecModal_wc, b0, env, fs, fVecModal_wc, 'fbfm');

if plotSpectrograms == 1
    figure
    spectrogram(real(ySAPFModal1), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using theoretical steel plate modal frequencies spectrogram 1')

    figure
    spectrogram(real(ySAPFModal2), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using theoretical steel plate modal frequencies spectrogram 2')

    figure
    spectrogram(real(ySAPFModal3), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using theoretical steel plate modal frequencies spectrogram 3')

    figure
    spectrogram(real(ySAPFModal4), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using theoretical steel plate modal frequencies spectrogram 4') 
end

%% TIME-VARYING APF vs. MODAL SYNTHESIS

% EXAMPLE 1: time-varying APF center/sounding frequencies = mesh modal 
% frequencies, fixed APF parameters
[yTVAPFModal1, yTVAPFModalMat1, TVAPFParams1] = TVAPFSynthesis(fVecModal, env, TVAPFParams, 0, fs);

% EXAMPLE 2: time-varying APF center/sounding frequencies = mesh modal 
% frequencies, randomized APF parameters - this sounds different everytime
[yTVAPFModal2, yTVAPFModalMat2, TVAPFParams2] = TVAPFSynthesis(fVecModal, env, TVAPFParams, 1, fs);

if plotSpectrograms == 1
    figure
    spectrogram(real(yTVAPFModal1), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying APF synthesis using theoretical steel plate modal frequencies spectrogram - fixed params')

    figure
    spectrogram(real(yTVAPFModal2), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying APF synthesis using theoretical steel plate modal frequencies spectrogram - random params')
end


%% STRETCHED APF & TIME-VARYING APF vs. MODAL SYNTHESIS
% all time-varying APF calls here use randomized parameters

b0 = (sqrt(1-B^2) - 1)/B;
fVecModal_wc = fVecModal*sqrt(1 - B^2);

% EXAMPLE 1: wc = mesh modal frequencies
[ySTVAPFModal1, ySTVAPFModalMat1, ~] = stretchedAPFAndTVAPFSynthesis(fVecModal, b0, env, TVAPFParams, 1, fs, [], 'none');

% EXAMPLE 2: w0 = mesh modal frequencies
[ySTVAPFModal2, ySTVAPFModalMat2, ~] = stretchedAPFAndTVAPFSynthesis(fVecModal_wc, b0, env, TVAPFParams, 1, fs, [], 'none');

% EXAMPLE 3: feedback FM with FBFM pitch glide with wc = mesh modal frequencies
[ySTVAPFModal3, ySTVAPFModalMat3, ~] = stretchedAPFAndTVAPFSynthesis(fVecModal, b0, env, TVAPFParams, 1, fs, fVecModal, 'fbfm');

% EXAMPLE 4: feedback FM with pitch glide with w0 = mesh modal frequencies
[ySTVAPFModal4, ySTVAPFModalMat4, ~] = stretchedAPFAndTVAPFSynthesis(fVecModal, b0, env, TVAPFParams, 1, fs, fVecModal_wc, 'fbfm');

if plotSpectrograms == 1
    figure
    spectrogram(real(ySTVAPFModal1), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using theoretical steel plate modal frequencies spectrogram 1')

    figure
    spectrogram(real(ySTVAPFModal2), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using theoretical steel plate modal frequencies spectrogram 2')

    figure
    spectrogram(real(ySTVAPFModal3), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using theoretical steel plate modal frequencies spectrogram 3')

    figure
    spectrogram(real(ySTVAPFModal4), hann(256), 128, 1024, fs, 'yaxis');
    title('stretched APF synthesis using theoretical steel plate modal frequencies spectrogram 4') 
end


%% write sounds to disk

if writeAudioFiles == 1
    
    outputDir = 'audioExamples/';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir)
    end
    
    % mesh
    audiowrite([outputDir 'yMeshReflCoeff.wav'], scaleForSavingAudio(yMeshReflCoeff), fs)
    audiowrite([outputDir 'yMeshFilter.wav'], scaleForSavingAudio(yMeshFilter), fs)
    
    % feedback FM vs. mesh
    audiowrite([outputDir 'yFBFMMesh1.wav'], scaleForSavingAudio(real(yFBFMMesh1)), fs)
    audiowrite([outputDir 'yFBFMMesh2.wav'], scaleForSavingAudio(real(yFBFMMesh2)), fs)
    audiowrite([outputDir 'yFBFMMesh3.wav'], scaleForSavingAudio(real(yFBFMMesh3)), fs)
    audiowrite([outputDir 'yFBFMMesh4.wav'], scaleForSavingAudio(real(yFBFMMesh4)), fs)
    
    % stretched APF vs. mesh
    audiowrite([outputDir 'ySAPFMesh1.wav'], scaleForSavingAudio(real(ySAPFMesh1)), fs)
    audiowrite([outputDir 'ySAPFMesh2.wav'], scaleForSavingAudio(real(ySAPFMesh2)), fs)
    audiowrite([outputDir 'ySAPFMesh3.wav'], scaleForSavingAudio(real(ySAPFMesh3)), fs)
    audiowrite([outputDir 'ySAPFMesh4.wav'], scaleForSavingAudio(real(ySAPFMesh4)), fs)
    
    % time-varying APF vs. mesh
    audiowrite([outputDir 'yTVAPFMesh1.wav'], scaleForSavingAudio(real(yTVAPFMesh1)), fs)
    audiowrite([outputDir 'yTVAPFMesh2.wav'], scaleForSavingAudio(real(yTVAPFMesh2)), fs)
    
    % stretched APF and time-varying APF vs. mes
    audiowrite([outputDir 'ySTVAPFMesh1.wav'], scaleForSavingAudio(real(ySTVAPFMesh1)), fs)
    audiowrite([outputDir 'ySTVAPFMesh2.wav'], scaleForSavingAudio(real(ySTVAPFMesh2)), fs)
    audiowrite([outputDir 'ySTVAPFMesh3.wav'], scaleForSavingAudio(real(ySTVAPFMesh3)), fs)
    audiowrite([outputDir 'ySTVAPFMesh4.wav'], scaleForSavingAudio(real(ySTVAPFMesh4)), fs)
    
    % steel plate modal synthesis
    audiowrite([outputDir 'yMS.wav'], scaleForSavingAudio(real(yMS)), fs)
    
    % feedback FM vs. steel plate modal synthesis
    audiowrite([outputDir 'yFBFMModal1.wav'], scaleForSavingAudio(real(yFBFMModal1)), fs)
    audiowrite([outputDir 'yFBFMModal2.wav'], scaleForSavingAudio(real(yFBFMModal2)), fs)
    audiowrite([outputDir 'yFBFMModal3.wav'], scaleForSavingAudio(real(yFBFMModal3)), fs)
    audiowrite([outputDir 'yFBFMModal4.wav'], scaleForSavingAudio(real(yFBFMModal4)), fs)
    
    % stretched APF vs. steel plate modal synthesis
    audiowrite([outputDir 'ySAPFModal1.wav'], scaleForSavingAudio(real(ySAPFModal1)), fs)
    audiowrite([outputDir 'ySAPFModal2.wav'], scaleForSavingAudio(real(ySAPFModal2)), fs)
    audiowrite([outputDir 'ySAPFModal3.wav'], scaleForSavingAudio(real(ySAPFModal3)), fs)
    audiowrite([outputDir 'ySAPFModal4.wav'], scaleForSavingAudio(real(ySAPFModal4)), fs)
    
    % time-varying APF vs. steel plate modal synthesis
    audiowrite([outputDir 'yTVAPFModal1.wav'], scaleForSavingAudio(real(yTVAPFModal1)), fs)
    audiowrite([outputDir 'yTVAPFModal2.wav'], scaleForSavingAudio(real(yTVAPFModal2)), fs)
    
    % stretched APF and time-varying APF vs. steel plate modal synthesis
    audiowrite([outputDir 'ySTVAPFModal1.wav'], scaleForSavingAudio(real(ySTVAPFModal1)), fs)
    audiowrite([outputDir 'ySTVAPFModal2.wav'], scaleForSavingAudio(real(ySTVAPFModal2)), fs)
    audiowrite([outputDir 'ySTVAPFModal3.wav'], scaleForSavingAudio(real(ySTVAPFModal3)), fs)
    audiowrite([outputDir 'ySTVAPFModal4.wav'], scaleForSavingAudio(real(ySTVAPFModal4)), fs)
    
    
    

end
