% This file creates figures for the SMC 2019 paper
outDir = '~/Documents/ucsd/winter2019/smc/background_planning/figure_making/';

fs = 44100;
dur = 1;

N = dur*fs;
T = 1/fs;
t = 0:T:(dur-T);


%% EXCITATION
%%%%%%%%%%%%%%%%%

attackTime = 0.07;
sustainTime = 0.03;
releaseTime = 0.15;

attackSamps = floor(attackTime * fs);
sustainSamps = floor(sustainTime * fs);
releaseSamps = floor(releaseTime * fs);

excitation = zeros(1, N);

excitation(1:attackSamps) = linspace(0, 1, attackSamps);
excitation(attackSamps+1:attackSamps+1+sustainSamps) = 1;
excitation(attackSamps+sustainSamps+1:attackSamps+sustainSamps+releaseSamps) = linspace(1, 0, releaseSamps);

figure
plot(t, excitation, 'linewidth', 2);
set(gca,'linewidth', 3)
set(gca,'XTick',[], 'YTick', [])
fig = gcf
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 1.6];
print([outDir 'excitation'], '-deps', '-r0')

%% EXCITATION DERIVATIVE
%%%%%%%%%%%%%%%%%%%%%%%%

dexcitation = [diff(excitation) 0];

figure
plot(t, dexcitation, 'linewidth', 2);
set(gca,'linewidth', 3)
set(gca,'XTick',[], 'YTick', [])
fig = gcf
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 1.6];
print([outDir 'excitation_derivative'], '-deps', '-r0')

%% LOOPBACK FM ADDITIVE
%%%%%%%%%%%%%%%%%%%%%%%%

% look at stretchedAPFSynthesis.m in proofOfConcept directory and print
% from there (envelopes are in there in my new update)

Nf = 3;
oscN = 500;
osc = zeros(Nf, oscN);
fVec = [233, 468, 706];
n = linspace(0, oscN-1, oscN);
b0 = -0.5;
for f=1:Nf
    osc(f,:) = (b0 + exp(j*2*pi*fVec(f)*n*T)) ./ (1 + b0*exp(j*2*pi*fVec(f)*n*T));
end


% LOOPBACK FM WITH ADDITIVE SYNTHESIS
figure
subplot(311); plot(real(osc(1, 1:oscN)), 'linewidth', 2)
set(gca,'linewidth', 3)
set(gca,'XTick',[], 'YTick', [])

subplot(312); plot(real(osc(2, 1:oscN)), 'linewidth', 2)
set(gca,'linewidth', 3)
set(gca,'XTick',[], 'YTick', [])

subplot(313); plot(real(osc(3, 1:oscN)), 'linewidth', 2)
set(gca,'linewidth', 3)
set(gca,'XTick',[], 'YTick', [])

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
print(['~/Documents/ucsd/winter2019/smc/background_planning/figure_making/' 'loopbackFMOsc2'], '-deps', '-r0')


%% RESONANT BODY IMPULSE RESPONSE
%%%%%%%%%%%%%%%%%
resIRWav = '../proofOfConcept/resonatorIRs/taiko/taiko2.wav';
[resIR, fsr] = audioread(resIRWav);
resIR = resIR(2159:end,1);
t = linspace(0, length(resIR)/fsr, length(resIR));

figure
plot(t, resIR, 'linewidth', 1);
xlim([t(1) t(end)])
set(gca,'linewidth', 3)
set(gca,'XTick',[], 'YTick', [])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
print([outDir 'resonantBodyIR'], '-deps', '-r0')


%% OUTPUT SIGNAL
%%%%%%%%%%%%%%%%%

addpath(genpath('../proofOfConcept'));

% additive synthesis x loopback FM signal
yMSWav = 'audioExamples/FBFMStretchedAPF/ySAPFModal1.wav';
[yMS, ~] = audioread(yMSWav);

e = zeros(N, 1);
attackSamps = 16;
sustainSamps = 2;
releaseSamps = 1024;
e(1:attackSamps) = linspace(0, 1, attackSamps);
e(attackSamps+1:attackSamps+1+sustainSamps) = 1;
e(attackSamps+sustainSamps+1:attackSamps+sustainSamps+releaseSamps) = linspace(1, 0, releaseSamps);

de = [diff(e); 0]; % take derivative for velocity?
y = percSynthExcitation(de, yMSWav);
y = y(1:N);

y = y/max(y);

figure
plot(t, y)
xlim([t(attackSamps+sustainSamps+releaseSamps) t(end)]);
set(gca,'linewidth', 3)
set(gca,'XTick',[], 'YTick', [])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
print([outDir 'percSynthesisOutput'], '-deps', '-r0')
