% plot decaying envelopes for SMC 2019 paper

%% input parameters 

% general
fs = 44100;
dur = 2;
plotSpectrograms = 0;

% feedback FM
B = 0.3;    % feedback coefficient
g = 0.9999; % pitch glide coefficient


%% derived parameters
N = fs*dur;
T = 1/fs;
BVec = g.^(0:N-1)';             % feedback FM pitch glide coefficient

%% ideal bar modal frequencies
% from Science of Percussion Instruments page 47

f1 = 440;
Nf = 7;
fVecBar = zeros(1, Nf);
for i=1:Nf
    if i==1
        fVecBar(i) = 3.011^2;
    else
        fVecBar(i) = (2*i+1)^2;
    end
end
fVecBar = fVecBar/fVecBar(1);
fVecBar = f1 * fVecBar;


% center frequencies for feedback FM if sounding frequencies = membrane
% modal freqeuncies
fcVecBar = fVecBar./(sqrt(1-B^2));

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

figure
for i=1:Nf
    plot(linspace(0, dur, N), env(i,:), 'linewidth', 2);
    hold on
end
xlabel('Time (sec)');
ylabel('Amplitude (linear)');
title('Decaying exponential amplitude envelopes');
set(gca, 'FontSize', 15)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
grid on
saveas(fig, 'figures/amplitudeEnvelopes', 'epsc')
