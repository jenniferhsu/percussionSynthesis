fs = 44100;
dur = 1;

T = 1/fs;
w = 2*pi*2000;

outDir = 'figures/';

% intial position
z = 1;
[hz,hp,ht] = zplane(z)
set(ht,'LineWidth', 2)
set(ht,'LineStyle', ':')
set(ht,'Color', [0.3 0.3 0.3])
set(hz, 'MarkerSize', 45)
set(hz, 'Marker', '.')
xlabel('')
ylabel('')
xlim([-1 1])
ylim([-1 1])
saveas(gcf, [outDir '/constantFreqOscillator0'], 'png')

% rotation
z = z*exp(j*w*T);
[hz,hp,ht] = zplane(z)
set(ht,'LineWidth', 2)
set(ht,'LineStyle', ':')
set(ht,'Color', [0.3 0.3 0.3])
set(hz, 'MarkerSize', 45)
set(hz, 'Marker', '.')
xlabel('')
ylabel('')
xlim([-1 1])
ylim([-1 1])
saveas(gcf, [outDir 'constantFreqOscillator1'], 'png')

% rotation
z = z*exp(j*w*T);
[hz,hp,ht] = zplane(z)
set(ht,'LineWidth', 2)
set(ht,'LineStyle', ':')
set(ht,'Color', [0.3 0.3 0.3])
set(hz, 'MarkerSize', 45)
set(hz, 'Marker', '.')
xlabel('')
ylabel('')
xlim([-1 1])
ylim([-1 1])
saveas(gcf, [outDir 'constantFreqOscillator2'], 'png')

% rotation
z = z*exp(j*w*T);
[hz,hp,ht] = zplane(z)
set(ht,'LineWidth', 2)
set(ht,'LineStyle', ':')
set(ht,'Color', [0.3 0.3 0.3])
set(hz, 'MarkerSize', 45)
set(hz, 'Marker', '.')
xlabel('')
ylabel('')
xlim([-1 1])
ylim([-1 1])
saveas(gcf, [outDir 'constantFreqOscillator3'], 'png')