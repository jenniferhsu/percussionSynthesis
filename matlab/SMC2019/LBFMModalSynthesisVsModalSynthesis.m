% LBFMModalSynthesisVsModalSynthesis.m
% This creates the figures that compare the results of Loopback FM Modal
% Synthesis with canonical modal synthesis (it would be additive synthesis
% in this case)

addpath(genpath('../modalFrequencyEquations'));
addpath(genpath('../proofOfConcept'));

%% input parameters 

% general
fs = 44100;
dur = 1;
plotSpectrograms = 0;
plotSMCFigures = 1;
saveSMCFigures = 0;

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

%% traditional modal/additive synthesis with exponential decay using the modes
yMS_low = zeros(1, N);
nVec = 0:T:(dur-T);

for i=1:Nf
    f = fVecModal_low(i);
    yMS_low = yMS_low + (exp(1j*2*pi*f*nVec) .* env);
end

if plotSpectrograms == 1
    figure
    spectrogram(real(yMS_low), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
end

%% Feedback FM

%%%%%%%%%%%%%%%%%
% low frequencies
%%%%%%%%%%%%%%%%%

% feedback FM sounding frequencies = modal frequencies
[yFBFMModal_low1, ~] = feedbackFMSynthesis(fcVecModal_low, B, env, fs);

% feedback FM with pitch glide with sounding frequencies = modal frequencies
[yFBFMModal_low2, ~] = feedbackFMSynthesis(fcVecModal_low, BVec, env, fs);

%%%%%%%%%%%%%%%%%%
% high frequencies
%%%%%%%%%%%%%%%%%%

% feedback FM sounding frequencies = modal frequencies
[yFBFMModal_high1, ~] = feedbackFMSynthesis(fcVecModal_high, B, env, fs);

% feedback FM with pitch glide with sounding frequencies = modal frequencies
[yFBFMModal_high2, ~] = feedbackFMSynthesis(fcVecModal_high, BVec, env, fs);


%% STRETCHED APF vs. MODAL SYNTHESIS

b0 = (sqrt(1-B^2) - 1)/B;

%%%%%%%%%%%%%%%%%
% low frequencies
%%%%%%%%%%%%%%%%%

% stretched APF, w0 = modal frequencies
%[ySAPFModal_low1, ~] = stretchedAPFSynthesis(fcVecModal_w0_low, b0, env, fs, [], 'none');
[ySAPFModal_low1, ~] = stretchedAPFSynthesis(fVecModal_low, b0, env, fs, [], 'none');

% stretched APF with pitch glide, w0 = modal frequencies
[ySAPFModal_low2, ~] = stretchedAPFSynthesis(fVecModal_low, b0, env, fs, fcVecModal_low, 'fbfm');

%%%%%%%%%%%%%%%%%%
% high frequencies
%%%%%%%%%%%%%%%%%%

% stretched APF, w0 = modal frequencies
[ySAPFModal_high1, ~] = stretchedAPFSynthesis(fcVecModal_high, b0, env, fs, [], 'none');

% stretched APF with pitch glide, w0 = modal frequencies
[ySAPFModal_high2, ~] = stretchedAPFSynthesis(fVecModal_high, b0, env, fs, fcVecModal_high, 'fbfm');


%% Musical Parameters - Sounding Frequency Diagram
% this figure compares how the spectrogram changes if we set the sounding
% frequencies to the modal frequencies vs. if we set the carrier frequencies
% to the modal frequencies

f_high2 = 5000;
fVecModal_high2 = modes*f_high2;
fcVecModal_high2 = fVecModal_high2./(sqrt(1-B^2));

[yFBFM_fcIsModal, ~] = feedbackFMSynthesis(fVecModal_high2, B, env, fs);

% feedback FM with pitch glide with sounding frequencies = modal frequencies
[yFBFM_f0IsModal, ~] = feedbackFMSynthesis(fcVecModal_high2, B, env, fs);


%% Figures for SMC

if plotSMCFigures==1
    
    % traditional modal synthesis vs. loopback FM modal synthesis with f_low=2000
    figure
%     subplot(311)
    spectrogram(real(yMS_low), hann(256), 128, 1024, fs, 'yaxis');
    title('Traditional modal synthesis spectrogram, f_c=2000Hz')
    set(gca,'FontSize',15)
    ylim([0 6]);
%     subplot(312)
%     spectrogram(real(yFBFMModal_low1), hann(256), 128, 1024, fs, 'yaxis');
%     title('feedback FM modal synthesis spectrogram, f_c=2000Hz')
%     set(gca,'FontSize',15)
%     subplot(313)
%     spectrogram(real(yFBFMModal_low2), hann(256), 128, 1024, fs, 'yaxis');
%     title('feedback FM modal synthesis with pitch glide spectrogram, f_c=2000Hz')
%     set(gca,'FontSize',15)
    if saveSMCFigures==1
        fig = gcf;
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 6 2.5];
        saveas(gcf, [saveDir 'traditionalMS'], 'epsc')
    end
    
    % loopback FM rotation vs stretched APF with f_low=2000 (static)
    figure
    subplot(211)
    spectrogram(real(yFBFMModal_low1), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM MS spectrogram, f_c=2000Hz')
    set(gca,'FontSize',15)
    subplot(212)
    spectrogram(real(ySAPFModal_low1), hann(256), 128, 1024, fs, 'yaxis');
    title('Stretched APF MS spectrogram, f_c=2000Hz')
    set(gca,'FontSize',15)
    if saveSMCFigures==1
        saveas(gcf, [saveDir 'rotationVsSAPF_f_low'], 'epsc')
    end
    
    % loopback FM rotation vs stretched APF with f_low=2000 (pitch glide)
    figure
    subplot(211)
    spectrogram(real(yFBFMModal_low2), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM MS with pitch glide spectrogram, f_c=2000Hz')
    set(gca,'FontSize',15)
    subplot(212)
    spectrogram(real(ySAPFModal_low2), hann(256), 128, 1024, fs, 'yaxis');
    title('Stretched APF MS with pitch glide spectrogram, f_c=2000Hz')
    set(gca,'FontSize',15)
    if saveSMCFigures==1
        saveas(gcf, [saveDir 'rotationVsSAPF_pitchGlide_f_low'], 'epsc')
    end
    
    % loopback FM rotation vs stretched APF with f_high=4000Hz (static)
    figure
    subplot(211)
    spectrogram(real(yFBFMModal_high1), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM MS spectrogram, f_c=4000Hz')
    set(gca,'FontSize',15)
    subplot(212)
    spectrogram(real(ySAPFModal_high1), hann(256), 128, 1024, fs, 'yaxis');
    title('Stretched APF MS spectrogram, f_c=4000Hz')
    set(gca,'FontSize',15)
    if saveSMCFigures==1
        saveas(gcf, [saveDir 'rotationVsSAPF_f_high'], 'epsc')
    end
    
    % loopback FM rotation vs stretched APF with f_high=4000Hz (pitch glide)
    figure
    subplot(211)
    spectrogram(real(yFBFMModal_high2), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM MS with pitch glide spectrogram, f_c=4000Hz')
    set(gca,'FontSize',15)
    subplot(212)
    spectrogram(real(ySAPFModal_high2), hann(256), 128, 1024, fs, 'yaxis');
    title('Stretched APF MS with pitch glide spectrogram, f_c=4000Hz')
    set(gca,'FontSize',15)
    if saveSMCFigures==1
        saveas(gcf, [saveDir 'rotationVsSAPF_pitchGlide_f_high'], 'epsc')
    end
    
    % sounding frequency (carrier or sounding frequencies set to equal
    % modal frequencies)
    figure
    subplot(211)
    spectrogram(real(yFBFM_fcIsModal), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM modal synthesis with f_c = f(i)')
    set(gca,'FontSize',15)
    subplot(212)
    spectrogram(real(yFBFM_f0IsModal), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM modal synthesis with f_0 = f(i)')
    set(gca,'FontSize',15)
    if saveSMCFigures==1
        saveas(gcf, [saveDir '?carrierOrSoundingAsModalFreq'], 'epsc')
    end
    
end



    