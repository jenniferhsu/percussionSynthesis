% pitch glides are not continuous with FM

fs = 44100;
dur = 1;

T = 1/fs;
N = fs*dur;

fc = 600;
fm = 1000;
B = 1;

BVec = linspace(0, 20, N);

y = zeros(1, N);
y(1) = 1;

for n=2:N
    %y = cos(2*pi*fc*n*T + B*sin(2*pi*fm*n*T));
    y(n) = cos(2*pi*fc*n*T + BVec(n)*sin(2*pi*fm*n*T));
end

figure
spectrogram(real(y), hann(256), 128, 1024, fs, 'yaxis');