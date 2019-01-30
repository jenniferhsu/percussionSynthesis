% percSynthTests.m
% this script shows how to use the percSynth function with different
% impulses

%% input parameters
fs = 44100;
N = 2*fs;

inputDir = 'audioExamples/membrane/';
%inputDir = 'audioExamples/circularPlate/';
%inputDir = 'audioExamples/woodBlocks/';

outputDir = ['audioExamples/' inputDir(15:end) 'commutedSynth/snareShell/'];
if ~exist(outputDir)
    mkdir(outputDir)
end

% resonating body impulse response locations
%resIRwav = 'resonatorIRs/CarpenterCenter.wav';
resIRwav = 'resonatorIRs/snareShell.wav';
%resIRwav = 'resonatorIRs/3000CStreetGarageStairwell.wav';

% instrument impulse response location
files = dir(inputDir);

for f=1:length(files)
    
    yWav = files(f).name;
    
    % skip files that are not .wav files
    if length(yWav) < 4
        continue;
    end
    if ~strcmp(yWav(end-2:end), 'wav')
        continue
    end

    %% create excitations and excitation derivatives
    
    excitations = zeros(N,4);
    d_excitations = zeros(N,4);

    %% 1. impulse
    excitations(1,1) = 1;
    d_excitations(:,1) = [diff(excitations(:,1)); 0];


    %% 2. attack, sustain, release signal
    % making attackSamps shorter retains more of the higher frequencies
    % increasing attackSamps to 512 removes the higher frequencies and creates
    % a very bassy sound.  attackSamps=8 sounds cool!
    %attackSamps = 512;
    attackSamps = 2;
    sustainSamps = 20;
    releaseSamps = 50;
    excitations(1:attackSamps, 2) = linspace(0, 1, attackSamps);
    excitations(attackSamps+1:attackSamps+1+sustainSamps, 2) = 1;
    excitations(attackSamps+sustainSamps+1:attackSamps+sustainSamps+releaseSamps, 2) = linspace(1, 0, releaseSamps);
    d_excitations(:,2) = [diff(excitations(:,2)); 0];

    %% 4. noise burst
    durNB = .10;
    sampNB = ceil(durNB*fs);
    excitations(:,3) = [(2*rand(1, sampNB))-1 zeros(1, N-sampNB)];
    d_excitations(:,3) = [diff(excitations(:,3)); 0];

    %% 5. filtered noise burst
    [B, A] = butter(5, [400/(fs/2) 2000/(fs/2)], 'bandpass');
    %freqz(B,A)
    excitations(:,4) = filter(B, A, excitations(:,3));
    d_excitations(:,4) = [diff(excitations(:,4)); 0];


    %% perform commuted synthesis

    y1 = percSynth(d_excitations(:,1), [inputDir yWav], resIRwav);
    y1 = y1(1:N);

    y2 = percSynth(d_excitations(:,2), [inputDir yWav], resIRwav);
    y2 = y2(1:N);

    y3 = percSynth(d_excitations(:, 3), [inputDir yWav], resIRwav);
    y3 = y3(1:N);

    y4 = percSynth(d_excitations(:,4), [inputDir yWav], resIRwav);
    y4 = y4(1:N);
    
    %% write to file
    audiowrite([outputDir yWav(1:end-4) '_impulse.wav'], scaleForSavingAudio(y1), fs);
    audiowrite([outputDir yWav(1:end-4) '_ASR.wav'], scaleForSavingAudio(y2), fs);
    audiowrite([outputDir yWav(1:end-4) '_noiseBurst.wav'], scaleForSavingAudio(y3), fs);
    audiowrite([outputDir yWav(1:end-4) '_noiseBurstFiltered.wav'], scaleForSavingAudio(y4), fs);
    
end
