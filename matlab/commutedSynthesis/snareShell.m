
% snare shell resonating mode frequencies
% from Rossing's Science of Percussion Instruments  (page 29)
f = [182, 330, 278, 341, 403, 442, 512, 556, 619];
F = length(f);

% modal frequency decay rates (dB/s) (page 32)
% drum stand decay rates (what do the ...s mean?)
d = [60, 30, 30, 30, 30, 50, 30, 35, 65];

% additive synthesis
dur = 2;
fs = 44100;
N = fs*dur;
T = 1/fs;
nT = 0:T:(dur-T);

y = zeros(1, N);

for i=1:F
    
    f0 = f(i);
    x = sin(2*pi*f0*nT);
    
    m = -d(i);
    envdB = m*nT;
    env = 10.^(envdB/20);
    
    xe = x.*env;
    
    y = y+xe;
end

audiowrite('resonatorIRs/snareShell.wav', scaleForSavingAudio(y), fs);