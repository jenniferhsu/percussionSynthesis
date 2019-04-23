% code to get the take the 3x3 transfer function of the digital waveguide
% mesh into a sum of first-order and second-order sections

% makes plots of:
%   1. the poles on the unit circle
%   2. a spectrogram of the the mesh

saveDissertationFigures = 1;

a = -0.999;
fs = 44100;
dur = 1;

N = fs*dur;
x = zeros(1, N);
x(1) = 1;

Nfft = 2^nextpow2(N);
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

% 3x3 mesh transfer function
B = [0 2 0 a 0 0 0 -a^2];
A = [1 0 (1+a)/2 0 0 0 -(a + a^2)/2 0 -a^2];
y = filter(B, A, x);

%% bank of 1st order sections
[r, p, k] = residuez(B,A);

y1 = zeros(size(y));
for i=1:length(r)
    BB = r(i);
    AA = [1 -p(i)];
    y1 = y1 + real(filter(BB, AA, x));
end

%% bank of 2nd order sections

% sort poles by angle
[~, pInds] = sort(angle(p));
ps = p(pInds);
rs = r(pInds);

% get poles with angle from 0 to pi
totNP = length(ps);
ps = ps(totNP/2:end);
rs = rs(totNP/2:end);

% form the 2nd order sections
Nsos = length(ps)-1;
sos = zeros(Nsos,6);
for i=1:Nsos
    
    if i==1
        
        % put the 0 and pi one into the sos first
        p1 = ps(1);
        p2 = ps(end);
        r1 = rs(1);
        r2 = rs(end);
        
        b0 = r1+r2; 
        b1 = -(r1*p2 + r2*p1);
        b2 = 0;
        a0 = 1;
        a1 = -(p1+p2);
        a2 = p1*p2;
    else
        
        b0 = 2*real(rs(i));
        b1 = -2*real(rs(i)*conj(ps(i)));
        b2 = 0;
        a0 = 1;
        a1 = -2*real(ps(i));
        a2 = abs(ps(i))^2;
    end
    
    sos(i,:) = [b0 b1 b2 a0 a1 a2];
    
end

y2 = zeros(size(y));
for i=1:Nsos
    BB = sos(i,1:3);
    AA = sos(i,4:6);
    xx = filter(BB, AA, x);
    y2 = y2 + xx;
    
    XX = fft(xx, Nfft);
    XXPos = XX(1:Nfft/2+1);
    
    plot(faxis, abs(XXPos), 'linewidth', 2);
    %xlabel('Frequency (Hz)');
    %ylabel('Amplitude (linear)');
    title(sprintf('p_%d=%.2f, r_%d=%.2f', i, ps(i), i, rs(i)))
    set(gca,'FontSize',15)
    xlim([faxis(1) faxis(end)]);
    ylim([0 4000]);
    grid on;
    if saveDissertationFigures==1
        fig = gcf;
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 6 2.5];
        saveas(gcf, ['figures/DWM3x3Biquad' num2str(i)], 'epsc')
    end
    
    %freqz(BB, AA);
    %pause
end

% Note that it can also be expressed as a series (product) of second order sections
% but that doesn't relate to how we're using it as modal synthesis

%% FFT check
Y = fft(y, Nfft);
YPos = Y(1:Nfft/2+1);

Y1 = fft(y1, Nfft);
Y1Pos = Y1(1:Nfft/2+1);

Y2 = fft(y2, Nfft);
Y2Pos = Y2(1:Nfft/2+1);


figure
subplot(311)
plot(y); grid on;
subplot(312)
plot(y1); grid on;
subplot(313)
plot(y2); grid on;
sgtitle('time-domain signal for 3x3 mesh filter and parallel filter banks');

figure
subplot(311)
plot(faxis, abs(YPos)); grid on;
subplot(312)
plot(faxis, abs(Y1Pos)); grid on;
subplot(313)
plot(faxis, abs(Y2Pos)); grid on;
sgtitle('magnitude spectrum of 3x3, center hit mesh and parallel filter bank');

%% pole-zero plot

figure
zplane([], p);
h = gca;
set(h.Children(2), 'MarkerSize', 15)
set(h.Children(2), 'LineWidth', 2)
set(h.Children(1), 'color', 'k')
xlabel('Real');
ylabel('Imaginary');
title('Poles for the 3x3 DWM on the unit circle');
set(gca,'FontSize',15)
if saveDissertationFigures==1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 8];
    print('figures/DWM3x3Poles', '-depsc', '-r0')
end

figure
plot(faxis, abs(YPos), 'linewidth', 2); 
xlim([0 faxis(end)])
set(gca,'FontSize',15)
xlabel('Frequency (Hz)');
ylabel('Amplitude (linear)');
title('3x3 DWM Magnitude Spectrum');
grid on;
if saveDissertationFigures==1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 2.5];
    saveas(gcf, ['figures/DWM3x3MagSpec'], 'epsc')
end
