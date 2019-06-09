fs = 44100;
dur = 1;

T = 1/fs;
N = fs*dur;

f0 = 600;
B = 0.5;

BVec = linspace(0, 0.999, N);

y = zeros(1, N);
y(1) = 1;
for n=2:N
    %y(n) = cos(2*pi*f0*n*T)*(1 + B*y(n-1));
    y(n) = cos(2*pi*f0*n*T)*(1 + BVec(n)*y(n-1));
end

Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

Y = fft(y, Nfft);
YPos = Y(1:Nfft/2+1);
figure
plot(faxis, abs(YPos))