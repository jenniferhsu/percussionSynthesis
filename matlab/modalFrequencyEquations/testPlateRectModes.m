% This script is a test file for plateRectModes(). It calls
% plateRect Modes and synthesizes a sound using the calculated
% modal frequencies.

%% inputs
dx = 51;     % plate size (x)
dy = 51;     % plate size (y)
h = 1;       % plate thickness
cL = 1;      % longitudinal wave speed
f0 = 440;    % fundamental frequency (Hz)
dur = 2;     % duration of hit in seconds
fs = 44100;  % sample rate (cycles/sec)


%% get the modal frequncies for a 51x51 sized plate
f = plateRectModes(dx, dy, h, cL);


%% synthesize a sound with the frequencies

% get the exact frequencies (Hz)
f = unique(f); 
f = f0*f;      	 

% exponentially decaying envelope
A = 1;			% initial amplitude
An = 0.00000001;	% final amplitude
l = -log(An)/(dur*fs);	% time constant

env = zeros(dur*fs, 1);
for t=1:(dur*fs)
    env(t) = A*exp(-l*t);
end

% synthesize sinusoids at the modal frequencies with the
% exponentially decaying envelope
y = zeros(dur*fs, 1);
tt = 0:1/fs:(dur-(1/fs));
mult = 1.0/length(f); % so we don't create a giant number
for i = 1:length(f)
    s = sin(2*pi*f0*f(i)*tt);
    y = y + mult*(env.*s');
end

soundsc(y,fs)