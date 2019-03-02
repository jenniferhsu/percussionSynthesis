addpath(genpath('../proofOfConcept/'))

fs = 44100;
dur = 1;

N = fs*dur;
T = 1/fs;
nT = 0:T:(dur-T);

f_c = 440;
x = sin(2*pi*f_c*nT);

% TVAPF parameters
f_pi = f_c; %   f_pi: the frequency at which the phase is -pi
f_m = 500;  %   f_m: frequency of modulation

% this is close to equal
f_b = 1;    %   f_b: the bandwidth of the generated sidebands
M = 1;      %   M: modulation depth
I = 0.01;   %   FM: index of modulation

f_b = 80;    %   f_b: the bandwidth of the generated sidebands
M = 700;      %   M: modulation depth
I = 1;   %   FM: index of modulation


% FFT parameters
Nfft = N;
faxis = (fs/2)*linspace(0, 1, Nfft/2+1);

%% FM

yfm = cos(2*pi*f_c*nT + I*sin(2*pi*f_m*nT));
YFM = fft(yfm, Nfft);
YFMPos = YFM(1:Nfft/2+1);
YFMPosAbsdB0 = 20*log10(abs(YFMPos)/max(abs(YFMPos)));

figure
plot(faxis, YFMPosAbsdB0, 'linewidth', 2)
hold on
plot([f_c f_c], [-60 0], 'r--');
% for i=1:5
%     fout1 = f_c+(i*f_m);
%     % aliased sidebands
%     if fout1 > fs/2
%         fout1 = (fs/2) - (fout - fs/2);
%     end
%     plot([fout1 fout1], [-60 0], 'g--');
%     
%     fout2 = f_c-(i*f_m);
%     % aliased sidebands
%     if fout2 < 0
%         fout2 = -1*fout2;
%     end
%     plot([fout2 fout2], [-60 0], 'm--');
% end
xlim([faxis(1) faxis(end)])
ylim([-60 0])
title('FM');

%% TVAPF
y = TVAPF2(x, f_pi, f_b, M, f_m, fs);
Y = fft(y, Nfft);
YPos = Y(1:Nfft/2+1);
YPosAbsdB = 20*log10(abs(YPos)/max(abs(YPos)));

figure
plot(faxis, YPosAbsdB, 'linewidth', 2)
hold on
plot([f_c f_c], [-60 0], 'r--');
% for i=1:5
%     fout1 = f_c+(i*f_m);
%     % aliased sidebands
%     if fout1 > fs/2
%         fout1 = (fs/2) - (fout - fs/2);
%     end
%     plot([fout1 fout1], [-60 0], 'g--');
%     
%     fout2 = f_c-(i*f_m);
%     % aliased sidebands
%     if fout2 < 0
%         fout2 = -1*fout2;
%     end
%     plot([fout2 fout2], [-60 0], 'm--');
% end
xlim([faxis(1) faxis(end)])
ylim([-60 0])
title('TV APF');

