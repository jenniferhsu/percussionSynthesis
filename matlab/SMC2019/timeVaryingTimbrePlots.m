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
g = 0.9999; % pitch/timbre glide coefficient

%% derived parameters
N = fs*dur;
T = 1/fs;
env = g.^(linspace(0, N, N));   % exponential decay
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient
b0Vec = (sqrt(1-BVec.^2) - 1) ./ BVec;

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
fcVecModal_low = fVecModal_low./(sqrt(1-BVec(end)^2));

%% Time-varying timbre for LBFM and SAPF

[yFBFMModal_tv, ~] = feedbackFMSynthesis(fcVecModal_low, BVec, env, fs);

[ySAPFModal_tv, ~] = stretchedAPFSynthesis(fVecModal_low, b0Vec, env, fs, [], 'none');

%% static timbre SAPF

[ySAPFModal_static, ~] = stretchedAPFSynthesis(fVecModal_low, b0Vec(1024), env, fs, [], 'none');


%% Figures for SMC

if plotSMCFigures==1
    
    % static vs time-varying timbre for SAPF
    figure
    subplot(311)
    spectrogram(real(ySAPFModal_static), hann(256), 128, 1024, fs, 'yaxis');
    title('Stretched APF MS, static timbre, b_0=-0.6312')
    set(gca,'FontSize',15)
    subplot(312)
    spectrogram(real(ySAPFModal_tv), hann(256), 128, 1024, fs, 'yaxis');
    title('Stretched APF MS, time-varying timbre, b(n)')
    set(gca,'FontSize',15)
    subplot(313)
    spectrogram(real(yFBFMModal_tv), hann(256), 128, 1024, fs, 'yaxis');
    title('Loopback FM MS, time-varying timbre, B(n)')
    set(gca,'FontSize',15)
    if saveSMCFigures==1
        fig = gcf;
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 6 8];
        print('figures/timeVaryingTimbre', '-depsc', '-r0')
    end
    
end



    