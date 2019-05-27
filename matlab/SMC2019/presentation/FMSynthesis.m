% FM synthesis example

addpath(genpath('../../proofOfConcept'))

% input parameters
fs = 44100;
dur = 1;
fc = 220;
fm = 1000;

% directories for saving
figDir = 'figures/';
audioDir = 'audioExamples/';

MVec = [0 1 5];
NM = length(MVec);



for m=1:NM
    M = MVec(m);

    % derived parameters
    T = 1/fs;
    N = dur*fs;
    n = 0:1:N-1;

    % fft parameters
    Nfft = 2^nextpow2(N);
    faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

    %% envelope
    fadeTime = 0.1; % seconds
    fadeSamps = fadeTime*fs;
    env = [ones(1, N-fadeSamps) linspace(1, 0, fadeSamps)];


    %% FM synthesis
    x = cos(2*pi*fc*n*T + M*sin(2*pi*fm*n*T));
    x = x .* env;

    X = fft(x, Nfft);
    XPos = X(1:Nfft/2+1);

    %% plots

    figure
    subplot(211)
    plot(linspace(0, dur, N), x, 'linewidth', 2);
    xlim([0, 0.01])
    xlabel('Time (seconds)');
    ylabel('Amplitude (linear)');
    set(gca, 'fontsize', 13);
    grid on

    subplot(212)
    plot(faxis, 20*log10(abs(XPos)/max(abs(XPos))), 'linewidth', 2)
    hold on
    plot([fc fc], [-60 0], 'r--', 'linewidth', 2)
    xlim([0 10000])
    ylim([-60 0])
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (dB)');
    set(gca, 'fontsize', 13);
    grid on
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 5];
    print([figDir 'FM_M' num2str(M) '.eps'], '-depsc', '-r0')
    
    %% save audio
    audiowrite([audioDir 'FM_M' num2str(M) '.wav'], scaleForSavingAudio(real(x)), fs);
    
    
end
