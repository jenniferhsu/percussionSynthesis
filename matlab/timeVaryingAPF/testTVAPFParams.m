fs = 44100;
dur = 1;

N = fs*dur;
T = 1/fs;
nT = 0:T:(dur-T);

f = 440;
x = sin(2*pi*f*nT);

% TVAPF parameters
f_pi = f;   %   f_pi: the frequency at which the phase is -pi
f_b = 100;  %   f_b: the bandwidth of the generated sidebands
M = 2000;   %   M: modulation depth
f_m = 500;  %   f_m: frequency of modulation

% FFT parameters
Nfft = N;
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

y = TVAPF2(x, f_pi, f_b, M, f_m, fs);
Y = fft(y, Nfft);
YPos = Y(1:Nfft/2+1);
YPosAbsdB = 20*log10(abs(YPos)/max(abs(YPos)));

figure
plot(faxis, YPosAbsdB, 'linewidth', 2)
hold on
plot([f f], [-60 0], 'r--');
for i=1:5
    fout1 = f+(i*f_m);
    % aliased sidebands
    if fout1 > fs/2
        fout1 = (fs/2) - (fout - fs/2);
    end
    plot([fout1 fout1], [-60 0], 'g--');
    
    fout2 = f-(i*f_m);
    % aliased sidebands
    if fout2 < 0
        fout2 = -1*fout2;
    end
    plot([fout2 fout2], [-60 0], 'm--');
end
xlim([faxis(1) faxis(end)])
ylim([-60 0])

figure
spectrogram(y, hann(256), 128, 1024, fs, 'yaxis');

%% What does f_b control?
% It looks like as f_b gets really large, the sideband levels begin to
% decrease. Setting f = 11000 and M=1000 makes it visible in the graphs, a
% setting like f=700 and M=500 is also pretty visible in the graphs.  With
% a setting like f=700, M=1000, it's much harder to summarize what is
% happening because there is a lot of aliasing, but you can hear a
% difference in the timbre. It sounds more nasally and then as the
% sidebands die away, it begins to sound muffled.
fbVec = [10, 100, 200, 500, 1000, 2000, 5000, 10000, 100000];
fbVec = [100:50:2000]
Nfb = length(fbVec);

yMat = zeros(N, Nfb);
YPosAbsdB = zeros(Nfft/2+1, Nfb);

f = 5500;
x = sin(2*pi*f*nT);
% x = zeros(N, 1);
% x(1) = 1;

% TVAPF parameters
f_pi = f;   %   f_pi: the frequency at which the phase is -pi
f_b = 100;  %   f_b: the bandwidth of the generated sidebands
M = 1000;    %   M: modulation depth
f_m = 500; %   f_m: frequency of modulation

figure
for i=1:Nfb
    yMat(:,i) = TVAPF2(x, f_pi, fbVec(i), M, f_m, fs);
    Y = fft(yMat(:,i), Nfft);
    YPos = Y(1:Nfft/2+1);
    YPosAbsdB(:,i) = 20*log10(abs(YPos)/max(abs(YPos)));
    
    % plot spectrogram
    
    hold on
    plot([f f], [-60 0], 'r--');    % frequency component at fHz
    for k=1:3                       % sidebands
        fout1 = f+(k*f_m);
        % aliased sidebands
        if fout1 > fs/2
            fout1 = (fs/2) - (fout - fs/2);
        end
        plot([fout1 fout1], [-100 0], 'g--');

        fout2 = f-(k*f_m);
        % aliased sidebands
        if fout2 < 0
            fout2 = -1*fout2;
        end
        plot([fout2 fout2], [-100 0], 'm--');
    end
    %plot(faxis, YPosAbsdB(:,i), 'linewidth', 2)
    plot(faxis, abs(YPos), 'linewidth', 2)
    hold off
    title(sprintf('%d', fbVec(i)))
    xlim([faxis(1) faxis(end)])
    ylim([0 24000]);
    %xlim([max(faxis(1),f-(2*f_m)) f+(2*f_m)])
    %ylim([-100 0])
    pause
    
end

%% How does M affect the signal

MVec = [100 1000 2000 4000];
NfM = length(MVec);

yMat = zeros(N, Nfb);
YPosAbsdB = zeros(Nfft/2+1, Nfb);

f = 700;%11000;  %700;
x = sin(2*pi*f*nT);
% x = zeros(N, 1);
% x(1) = 1;

% TVAPF parameters
f_pi = f;   %   f_pi: the frequency at which the phase is -pi
f_b = 100;  %   f_b: the bandwidth of the generated sidebands
M = 500;    %   M: modulation depth
f_m = 1000; %   f_m: frequency of modulation

figure
for i=1:NfM
    yMat(:,i) = TVAPF2(x, f_pi, f_b, MVec(i), f_m, fs);
    Y = fft(yMat(:,i), Nfft);
    YPos = Y(1:Nfft/2+1);
    YPosAbsdB(:,i) = 20*log10(abs(YPos)/max(abs(YPos)));
    
    % plot spectrogram
    
    hold on
    plot([f f], [-60 0], 'r--');    % frequency component at fHz
    for k=1:3                       % sidebands
        fout1 = f+(k*f_m);
        % aliased sidebands
        if fout1 > fs/2
            fout1 = (fs/2) - (fout - fs/2);
        end
        plot([fout1 fout1], [-100 0], 'g--');

        fout2 = f-(k*f_m);
        % aliased sidebands
        if fout2 < 0
            fout2 = -1*fout2;
        end
        plot([fout2 fout2], [-100 0], 'm--');
    end
    %plot(faxis, YPosAbsdB(:,i), 'linewidth', 2)
    plot(faxis, abs(YPos), 'linewidth', 2)
    hold off
    title(sprintf('%d', MVec(i)))
    xlim([faxis(1) faxis(end)])
    %ylim([-100 0])
    pause
    
end










