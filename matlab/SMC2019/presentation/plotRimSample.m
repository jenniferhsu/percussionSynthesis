[y,fs] = audioread('/Users/jenniferhsu/Documents/ucsd/dissertation/presentations/SMC2019/sound examples/rim-original.wav');

y = y(:,1);
figDir = 'figures/';


N = length(y);
dur = N/fs;
t = linspace(0, dur, N);

figure
plot(t, y, 'linewidth', 2);
xlabel('Time (seconds)');
ylabel('Amplitude (linear)');
xlim([0 dur]);
grid on
set(gca, 'fontsize', 15);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
print([figDir 'rim-original.eps'], '-depsc', '-r0')