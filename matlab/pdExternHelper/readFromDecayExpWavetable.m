% how to read from a decaying exponential table

% read in wavetable to a vector
fileID = fopen('decayExpWavetable.txt','r');
x = fscanf(fileID, '%f');
N = length(x);      % length of wavetable
T60EN = 0.75 * N;    % envelope wavetable T60;

% sample rate and signal settings
fs = 44100;
dur = 2;
T = 1/fs;

% desired T60 for envelope (in seconds)
T60 = 0.8;

%% use wavetable to create a sine wave
T60N = T60 * fs; % T60 in samples
einc = T60EN / T60N;

y = zeros(1, dur*fs);
ei = 0;
eind0 = 0;
eind1 = 0;

for i=1:dur*fs
    
    if ei >= N-1
        y(i) = 0;
    else
        % linear interpolation
        ind0 = floor(ei);
        frac = ei - ind0;

        % index by 1
        ind0 = ind0 + 1;
        ind1 = ind0 + 1;

        % perform linear interpolation
        y(i) = x(ind0) * (1 - frac) + x(ind1) * frac;

        ei = ei + einc;
    end
 
end

% compare with ysin to make sure it's correct
A = 1;
tau = -(T60 * fs)/log(0.001);
yenv = A * exp(-(0:(fs*dur-1))/tau);

figure
plot(y)
hold on
plot(yenv)

sum(y - yenv)

figure
plot(linspace(0, dur, dur*fs), 20*log10(y/max(y)));
hold on
plot(linspace(0, dur, dur*fs), 20*log10(yenv/max(yenv)));
plot([T60 T60], [-100, 0], 'r')
xlim([0 dur])
ylim([-100 0])
