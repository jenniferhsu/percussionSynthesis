% create a 2048-sample long exponential decay table
N = 2048;
n = 0:N-1;

% make the exponentially decaying envelope where -60dB occurs at 3/4 of
% the way to the end of the envelope
A = 1;
N60 = 0.75 * N;   % T60 as a sample
tau = -N60 / log(0.001);

x = A * exp(-n/tau);

% check that it is what we expect
figure
plot(20*log10(x/max(x)));
hold on
plot([N60 N60], [-90 0], 'r')
xlim([0 N-1])

% write to file
fileID = fopen('decayExpWavetable.txt', 'w');
for i=1:N
    fprintf(fileID, '%.12f ', x(i));
end
fclose(fileID);
