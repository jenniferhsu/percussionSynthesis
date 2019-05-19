% envMat_Tests.m: tests that envMat.m is creating the correct amplitude
% envelopes
%
% author: Jennifer Hsu
% date: Spring 2019

e0Vec = [1, 0.9, 0.8, 0.7 0.6]';
T60Vec = [0.9, 0.8, 0.7, 0.6, 0.5]';
dur = 1.0;
fs = 44100;
env = envMat(e0Vec, T60Vec, dur, fs);

Nf = size(env, 1);
N = size(env, 2);
nT = 0:1/fs:(dur-1/fs);

figure
plot([0 dur], [-60 -60], 'r--')
hold on
for f=1:Nf
    plot(nT, 20*log10(abs(env(f,:))/max(abs(env(f,:)))));
    hold on
    plot([T60Vec(f) T60Vec(f)], [-100 0], 'r--')
end
ylim([-100 0])
title('decaying exponential amplitude envelopes and T60 times')