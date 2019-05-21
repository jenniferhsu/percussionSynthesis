% TVAPFModalSynthesisVsModalSynthesis.m
% This creates the figures that compares:
%   a. Traditional Modal Synthesis
%   b. Time-Varying APF Traditional Modal Synthesis
%   c. Loopback FM Modal Synthesis
%   d. Time-Varying APF Loopback FM Modal Synthesis

addpath(genpath('../modalFrequencyEquations'));
addpath(genpath('../proofOfConcept'));

%% input parameters 

% general
fs = 44100;
dur = 1;
plotSpectrograms = 0;
plotFigures = 1;
saveFigures = 1;

saveDir = 'figures/';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% feedback FM
f_low = 2000;
f_high = 4000;
B = 0.9;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

%% derived parameters
N = fs*dur;
T = 1/fs;
env = g.^(linspace(0, N, N));   % exponential decay
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient
nVec = 0:T:(dur-T);

%% Modal synthesis frequencies
% generate modes for a 3x3 mode steel plate using the von Karman equations

modes = plateRectModes(3, 3, 1, 1);
modes = modes(:);
modes = unique(modes);
Nf = length(modes);

% f_low = 2000
fVecModal_low = modes*f_low;
% carrier frequencies for feedback FM if sounding frequencies = von Karman 
% modal frequencies
fcVecModal_low = fVecModal_low./(sqrt(1-B^2));

% f_high = 5500;
fVecModal_high = modes*f_high;
fcVecModal_high = fVecModal_high./(sqrt(1-B^2));


%% traditional modal synthesis with exponential decay using the modes
yMS_low = zeros(1, N);

for i=1:Nf
    f = fVecModal_low(i);
    yMS_low = yMS_low + (exp(1j*2*pi*f*nVec) .* env);
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yMS_low), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
end


%% time-varying allpass filtered traditional modal synthesis with 
% exponential decay using the modes

yMSTVAPF_low = zeros(1, N);

f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 1000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

for i=1:Nf
    f_pi = fVecModal_low(i);
    x1 = exp(1j*2*pi*f_pi*nVec);
    y1 = TVAPF2(x1, f_pi, f_b, M, f_m, fs);
    yMSTVAPF_low = yMSTVAPF_low + (y1 .* env);
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yMSTVAPF_low), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying allpass filtered modal synthesis spectrogram')
end


%% Feedback FM

%%%%%%%%%%%%%%%%%
% low frequencies
%%%%%%%%%%%%%%%%%

% feedback FM sounding frequencies = modal frequencies
[yFBFMModal_low1, yFBFMModalMat_low1] = feedbackFMSynthesis(fcVecModal_low, B, env, fs);

% feedback FM with pitch glide with sounding frequencies = modal frequencies
[yFBFMModal_low2, yFBFMModalMat_low2] = feedbackFMSynthesis(fcVecModal_low, BVec, env, fs);

%%%%%%%%%%%%%%%%%%
% high frequencies
%%%%%%%%%%%%%%%%%%

% feedback FM sounding frequencies = modal frequencies
[yFBFMModal_high1, yFBFMModalMat_high1] = feedbackFMSynthesis(fcVecModal_high, B, env, fs);

% feedback FM with pitch glide with sounding frequencies = modal frequencies
[yFBFMModal_high2, yFBFMModalMat_high2] = feedbackFMSynthesis(fcVecModal_high, BVec, env, fs);


%% time-varying allpass filtered loopback FM modal synthesis (low)
% (no pitch glide)

yLBTVAPF_low1 = zeros(1, N);

f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 2000;   %   M: modulation depth
f_m = 1000;  %   f_m: frequency of modulation

for i=1:Nf
    f_pi = fVecModal_low(i);
    x1 = yFBFMModalMat_low1(i,:);
    y1 = TVAPF2(x1, f_pi, f_b, M, f_m, fs);
    yLBTVAPF_low1 = yLBTVAPF_low1 + y1;
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yLBTVAPF_low1), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying allpass filtered loopback FM spectrogram - f_c = 2000Hz')
end

%% time-varying allpass filtered loopback FM modal synthesis (low)
% (with pitch glide)

yLBTVAPF_low2 = zeros(1, N);

f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 2000;   %   M: modulation depth
f_m = 1000;  %   f_m: frequency of modulation

for i=1:Nf
    f_pi = fVecModal_low(i);
    x1 = yFBFMModalMat_low2(i,:);
    y1 = TVAPF2(x1, f_pi, f_b, M, f_m, fs);
    yLBTVAPF_low2 = yLBTVAPF_low2 + y1;
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yLBTVAPF_low2), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying allpass filtered loopback FM (with pitch glide) spectrogram - f_c = 2000Hz')
end

%% time-varying allpass filtered loopback FM modal synthesis (high)
% (no pitch glide)

yLBTVAPF_high1 = zeros(1, N);

f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 1000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

for i=1:Nf
    f_pi = fVecModal_high(i);
    x1 = yFBFMModalMat_high1(i,:);
    y1 = TVAPF2(x1, f_pi, f_b, M, f_m, fs);
    yLBTVAPF_high1 = yLBTVAPF_high1 + y1;
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yLBTVAPF_high1), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying allpass filtered loopback FM spectrogram (f_c = 4000Hz)')
end

%% time-varying allpass filtered loopback FM modal synthesis (low)
% (with pitch glide)

yLBTVAPF_high2 = zeros(1, N);

f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 1000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

for i=1:Nf
    f_pi = fVecModal_high(i);
    x1 = yFBFMModalMat_high2(i,:);
    y1 = TVAPF2(x1, f_pi, f_b, M, f_m, fs);
    yLBTVAPF_high2 = yLBTVAPF_high2 + y1;
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yLBTVAPF_high2), hann(256), 128, 1024, fs, 'yaxis');
    title('time-varying allpass filtered loopback FM (with pitch glide) spectrogram (f_c = 4000Hz)')
end



%% Figures

if plotFigures==1
    
    % time-varying allpass filter using sinusoids of traditional 
    % modal synthesis with f_low=2000
    figure
    spectrogram(real(yMSTVAPF_low), hann(256), 128, 1024, fs, 'yaxis');
    title('Spectrogram of time-varying allpass filtered MS, f_c=2000Hz')
    set(gca,'FontSize',15)
    ylim([0 8]);

    if saveFigures==1
        fig = gcf;
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 6 2.5];
        saveas(gcf, [saveDir 'timeVaryingAPFModalSynthesis'], 'epsc')
    end
    
    % loopback FM rotation vs TV APF of loopback FM with f_low=2000 (static)
    figure
    subplot(211)
    spectrogram(real(yFBFMModal_low1), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM MS spectrogram')
    set(gca,'FontSize',15)
    subplot(212)
    spectrogram(real(yLBTVAPF_low1), hann(256), 128, 1024, fs, 'yaxis');
    title('Time-varying APF of loopback FM MS spectrogram')
    set(gca,'FontSize',15)
    if saveFigures==1
        saveas(gcf, [saveDir 'loopbackFMVsTimeVaryingAPF_f_low'], 'epsc')
    end
    
    % loopback FM rotation vs TV APF of loopback FM with f_low=2000 (pitch glide)
    figure
    subplot(211)
    spectrogram(real(yFBFMModal_low2), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM MS with pitch glide spectrogram')
    set(gca,'FontSize',15)
    subplot(212)
    spectrogram(real(yLBTVAPF_low2), hann(256), 128, 1024, fs, 'yaxis');
    title('TV APF of loopback FM MS with pitch glide spectrogram')
    set(gca,'FontSize',15)
    if saveFigures==1
        saveas(gcf, [saveDir 'loopbackFMVsTimeVaryingAPF_pitchGlide_f_low'], 'epsc')
    end
    
    % loopback FM rotation vs TV APF of loopback FM with f_high=4000 (static)
    figure
    subplot(211)
    spectrogram(real(yFBFMModal_high1), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM MS spectrogram')
    set(gca,'FontSize',15)
    subplot(212)
    spectrogram(real(yLBTVAPF_high1), hann(256), 128, 1024, fs, 'yaxis');
    title('Time-varying APF of loopback FM MS spectrogram')
    set(gca,'FontSize',15)
    if saveFigures==1
        saveas(gcf, [saveDir 'loopbackFMVsTimeVaryingAPF_f_high'], 'epsc')
    end
    
    % loopback FM rotation vs TV APF of loopback FM with f_high=4000 (pitch glide)
    figure
    subplot(211)
    spectrogram(real(yFBFMModal_high2), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM MS with pitch glide spectrogram')
    set(gca,'FontSize',15)
    subplot(212)
    spectrogram(real(yLBTVAPF_high2), hann(256), 128, 1024, fs, 'yaxis');
    title('TV APF of loopback FM MS with pitch glide spectrogram')
    set(gca,'FontSize',15)
    if saveFigures==1
        saveas(gcf, [saveDir 'loopbackFMVsTimeVaryingAPF_pitchGlide_f_high'], 'epsc')
    end
    
  
end

