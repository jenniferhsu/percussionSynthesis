% closed-form loopback FM examples

addpath(genpath('../proofOfConcept'))

% input parameters
fs = 44100;
dur = 1;

% directories for saving
figDir = 'figures/';
audioDir = 'audioExamples/';

% derived parameters
N = fs*dur;
n = 0:1:N-1;
T = 1/fs;

% fft parameters
Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

col = colormap('pink');


%% envelope
fadeTime = 0.1; % seconds
fadeSamps = fadeTime*fs;
env = [ones(1, N-fadeSamps) linspace(1, 0, fadeSamps)];


%% changing w0

%f0Vec = [220, 261.63, 329.63, 392, 493.88, 587.33];
f0Vec = [220, 329.63, 493.88];
b0 = -0.5;
Nf = length(f0Vec);

z0Matf0 = zeros(Nf, N);
Z0Matf0Pos = zeros(Nf, Nfft/2+1);
for i=1:Nf
    % synthesize the signal
    w0 = 2*pi*f0Vec(i);
    z0Matf0(i,:) = (b0 + exp(1j*w0*n*T)) ./ (1 + b0*exp(1j*w0*n*T));
    
    z0Matf0(i,:) = z0Matf0(i,:) .* env;
    
    % fft
    Z = fft(real(z0Matf0(i,:)), Nfft);
    Z0Matf0Pos(i,:) = Z(1:Nfft/2+1);
end

%% changing b0

b0Vec = -[0 0.4 0.8];
f0 = 392;
w0 = 2*pi*f0;
Nb0 = length(b0Vec);

z0Matb0 = zeros(Nb0, N);
Z0Matb0Pos = zeros(Nb0, Nfft/2+1);
for i=1:Nb0
    b0 = b0Vec(i);
    z0Matb0(i,:) = (b0 + exp(1j*w0*n*T)) ./ (1 + b0*exp(1j*w0*n*T));
    
    z0Matb0(i,:) = z0Matb0(i,:) .* env;
    
    % fft
    Z = fft(real(z0Matb0(i,:)), Nfft);
    Z0Matb0Pos(i,:) = Z(1:Nfft/2+1);
end


%% plots

% changing w0 for static equation
for i=1:Nf
    f0 = f0Vec(i);
    
    figure
    subplot(211)
    plot(n*T, real(z0Matf0(i,:)), 'linewidth', 2);
    %xlabel('Time (seconds)');
    %ylabel('Amplitude (linear)');
    title('Time-domain waveform');
    xlim([0 0.1])
    set(gca, 'fontsize', 13);
    grid on
    
    subplot(212)
    plot(faxis, 20*log10(abs(Z0Matf0Pos(i,:))/max(abs(Z0Matf0Pos(i,:)))), 'linewidth', 2);
    %plot(faxis, abs(Z0Matf0Pos(i,:)), 'linewidth', 2);
    hold on
    plot([f0 f0], [-60 0], 'color', col(32,:), 'linestyle', '--', 'linewidth', 2);
    %xlabel('Frequency (Hz)')
    %ylabel('Magnitude (dB)')
    title('Magnitude spectrum');
    xlim([0 6000]);
    ylim([-60 0]);
    grid on
    set(gca, 'fontsize', 13);
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    print([figDir 'closedForm_f0' num2str(f0) '.eps'], '-depsc', '-r0')
    %close
end

% changing b0 for static equations
for i=1:Nb0
    b0 = b0Vec(i);
    
    figure
    subplot(211)
    plot(n*T, real(z0Matb0(i,:)), 'linewidth', 2);
    %xlabel('Time (seconds)');
    %ylabel('Amplitude (linear)');
    title('Time-domain waveform');
    xlim([0 0.01])
    set(gca, 'fontsize', 13);
    grid on
    
    subplot(212)
    plot(faxis, 20*log10(abs(Z0Matb0Pos(i,:))/max(abs(Z0Matb0Pos(i,:)))), 'linewidth', 2);
    hold on
    plot([392 392], [-60 0], 'color', col(32,:), 'linestyle', '--', 'linewidth', 2);
    %xlabel('Frequency (Hz)')
    %ylabel('Magnitude (dB)')
    title('Magnitude spectrum');
    xlim([0 12000]);
    ylim([-60 0]);
    grid on
    set(gca, 'fontsize', 13);
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    print([figDir 'closedForm_b0' num2str(b0) '.eps'], '-depsc', '-r0')
    close
end






%% save audio

% chanigng w0
for i=1:Nf
    f0 = f0Vec(i);
    audiowrite([audioDir 'closedForm_b0Static_f0' num2str(f0) '.wav'], scaleForSavingAudio(real(z0Matf0(i,:))), fs);
end

% changing b0
for i=1:Nf
    b0 = b0Vec(i);
    audiowrite([audioDir 'closedForm_f0Static_b0' num2str(b0) '.wav'], scaleForSavingAudio(real(z0Matb0(i,:))), fs);
end
