% how to read from a wavetable

% read in wavetable to a vector
fileID = fopen('sinWavetable.txt','r');
x = fscanf(fileID, '%f');
N = length(x); % length of wavetable

% sample rate and signal settings
fs = 44100;
dur = 1;
T = 1/fs;

% desired frequency for sine wave
f = 220;

%% use wavetable to create a sine wave
y = zeros(1, dur*fs);
p = zeros(1, dur*fs);
for i=1:dur*fs
    if i==1
        p(i) = 0;
    else
        p(i) = p(i-1) + (N*T*f);
    end
    p(i) = mod(p(i), 2048);
    
    % linear interpolation
    ind0 = floor(p(i));% + 1;
    frac = p(i) - ind0;
    
    % index by 1
    ind0 = ind0 + 1;
    if ind0 > N
        ind0 = 1;
    end
    ind1 = ind0+1;
    if ind1 > N
        ind1 = 1;
    end
    
    % perform the linear interpolation
    y(i) = x(ind0) * (1 - frac) + x(ind1) * frac;
    
end

% compare with ysin to make sure it's correct
ysin = sin(2*pi*f*(0:T:(1-T)));

%% use wavetable to create a cosine wave
y2 = zeros(1, dur*fs);
p = zeros(1, dur*fs);

for i=1:dur*fs
    if i==1
        p(i) = (N/4) - 0;
    else
        p(i) = p(i-1) + (N*T*f);
    end
    p(i) = mod(p(i), 2048);
    
    % linear interpolation
    ind0 = floor(p(i));% + 1;
    frac = p(i) - ind0;
    
    % index by 1
    ind0 = ind0 + 1;
    if ind0 > N
        ind0 = 1;
    end
    ind1 = ind0+1;
    if ind1 > N
        ind1 = 1;
    end
    
    % perform the linear interpolation
    y2(i) = x(ind0) * (1 - frac) + x(ind1) * frac;
    
end

% compare with ycos to make sure it's correct
ycos = cos(2*pi*f*(0:T:(1-T)));
