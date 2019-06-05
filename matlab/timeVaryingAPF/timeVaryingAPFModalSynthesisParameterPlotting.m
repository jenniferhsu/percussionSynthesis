%% timeVaryingAPFModalSynthesisParameterPlotting.m
% This script plots how changing the parameters in the time-varying allpass
% filter affects the modal synthesis results for traditional modal
% synthesis with sinusoids

addpath(genpath('../modalFrequencyEquations'));
addpath(genpath('../proofOfConcept'));

%% input parameters 

% general
fs = 44100;
dur = 1;
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

if plotFigures == 1
    figure
    spectrogram(real(yMS_low), hann(256), 128, 1024, fs, 'yaxis');
    title('modal synthesis spectrogram')
end

%% time-varying allpass filtered traditional modal synthesis with 
% exponential decay using the modes
% THE EFFECT OF CHANGING M

M = 1000;   %   M: modulation depth
MVec = [100 500 1000 2000 5000 10000 20000 50000];
NM = length(MVec);

yMSTVAPF_M = zeros(NM, N);

f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
f_m = 500;  %   f_m: frequency of modulation

for m=1:NM
    M = MVec(m);
    for i=1:Nf
        f_pi = fVecModal_low(i);
        x1 = exp(1j*2*pi*f_pi*nVec);
        y1 = TVAPF2(x1, f_pi, f_b, M, f_m, fs);
        yMSTVAPF_M(m,:) = yMSTVAPF_M(m,:) + (y1 .* env);
    end
end

if plotFigures == 1
    fig = figure
    for m=1:NM
        subplot(4, 2, m)
        spectrogram(real(yMSTVAPF_M(m,:)), hann(256), 128, 1024, fs, 'yaxis');
        colorbar('off')
        title(sprintf('M_i=%1.f', MVec(m)), 'fontsize', 15)
        ylim([0 20])
    end
    sgtitle('spectrogram of time-varying AP_i(n) MS for different M_i values', 'fontsize', 15)
    if saveFigures == 1
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 6 6];
        print('figures/timeVaryingAPFVaryingM', '-depsc', '-r0')
    end
end

%% time-varying allpass filtered traditional modal synthesis with 
% exponential decay using the modes
% THE EFFECT OF CHANGING f_m

f_m = 500;  %   f_m: frequency of modulation
fmVec = [10 100 200 500 1000 2000 3000 5000];
Nfm = length(fmVec);

yMSTVAPF_f_m = zeros(NM, N);

f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
M = 1000;   %   M: modulation depth

for m=1:Nfm
    f_m = fmVec(m);
    for i=1:Nf
        f_pi = fVecModal_low(i);
        x1 = exp(1j*2*pi*f_pi*nVec);
        y1 = TVAPF2(x1, f_pi, f_b, M, f_m, fs);
        yMSTVAPF_f_m(m,:) = yMSTVAPF_f_m(m,:) + (y1 .* env);
    end
end

if plotFigures == 1
    figure
    for m=1:Nfm
        subplot(4, 2, m)
        spectrogram(real(yMSTVAPF_f_m(m,:)), hann(256), 128, 1024, fs, 'yaxis');
        colorbar('off')
        title(sprintf('f_{m,i}=%1.f', fmVec(m)), 'fontsize', 15)
        ylim([0 20])
    end
    sgtitle('spectrogram of time-varying AP_i(n) MS for different f_{m,i} values', 'fontsize', 15)
    if saveFigures == 1
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 6 6];
        print('figures/timeVaryingAPFVaryingfm', '-depsc', '-r0')
    end
end



%% time-varying allpass filtered traditional modal synthesis with 
% exponential decay using the modes
% THE EFFECT OF CHANGING f_b


f_b = 100;  %   f_b: the bandwidth of the generated sidebands?
fbVec = 10.^(0:5);
Nfb = length(fbVec);

yMSTVAPF = zeros(Nfb, N);

M = 1000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

for b=1:Nfb
    f_b = fbVec(b);
    for i=1:Nf
        f_pi = fVecModal_low(i);
        x1 = exp(1j*2*pi*f_pi*nVec);
        y1 = TVAPF2(x1, f_pi, f_b, M, f_m, fs);
        yMSTVAPF(b,:) = yMSTVAPF(b,:) + (y1 .* env);
    end
end

if plotFigures == 1
    figure
    for b=1:Nfb
        subplot(3, 2, b)
        spectrogram(real(yMSTVAPF(b,:)), hann(256), 128, 1024, fs, 'yaxis');
        colorbar('off')
        ylim([0 8])
        title(sprintf('f_{b,i}=%1.f', fbVec(b)))
    end
    sgtitle('spectrogram of time-varying AP_i(n) MS for different f_{b,i} values', 'fontsize', 15)
    if saveFigures == 1
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 6 6];
        print('figures/timeVaryingAPFVaryingfb', '-depsc', '-r0')
    end
end
