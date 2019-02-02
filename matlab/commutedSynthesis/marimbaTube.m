% closed tube resonating frequencies
% from Rossing's Science of Percussion Instruments  (page 50)

addpath(genpath('../proofOfConcept/'));

dur = 2;
fs = 44100;
N = fs*dur;
T = 1/fs;
nT = 0:T:(dur-T);

%% modal frequencies
Nf = 9;
f = zeros(1, Nf);
m = 1;
f(1) = 223;
for i=2:Nf
    m = m + 2;
    f(i) = m*f(1);
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

g = 0.9999;
e = g.^(linspace(0, N, N));   % exponential decay

env = zeros(Nf, N);
for i=1:Nf
    env(i,:) = aStart(i) * e;
end

%% additive synthesis

y = zeros(1, N);

for i=1:Nf
    
    f0 = f(i);
    x = sin(2*pi*f0*nT);
    
    xe = x .* env(i,:);
    
    y = y+xe;
end

audiowrite('../proofOfConcept/resonatorIRs/marimbaTube.wav', scaleForSavingAudio(y), fs);