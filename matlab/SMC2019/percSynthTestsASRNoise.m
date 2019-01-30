% percSynthTestsSMC2019.m
% this script shows how to use the percSynth function with different
% impulses

addpath('../proofOfConcept');

%% input parameters
fs = 44100;
N = 2*fs;

WRITETOFILE = 0; % set to 1 to save audio examples
outputDir = 'audioExamples/excitations/';

% modal synthesis signal
yMSWav = '../proofOfConcept/audioExamples/timeVaryingAPF/yTVAPFModal2.wav';
%yMSWav = '../proofOfConcept/audioExamples/FBFMStretchedAPF/ySAPFModal1.wav';

% input percussion signal for comparison
[yMS, ~] = audioread(yMSWav);

% get things ready for saving
if ~exist([outputDir 'ASR/'])
    mkdir([outputDir 'ASR/'])
end
if ~exist([outputDir 'noiseBurst/'])
    mkdir([outputDir 'noiseBurst/'])
end

tokens = strsplit(yMSWav, '/');
wavFileName = tokens{end};



%% ATTACK, SUSTAIN, RELEASE ENVELOPES
% making attackSamps shorter retains more of the higher frequencies
% increasing attackSamps to 512 removes the higher frequencies and creates
% a very bassy sound.  attackSamps=8 sounds cool!

attackSampsVec = [2, 20, 200, 2000];
sustainSampsVec = [2, 20, 200, 2000];
releaseSampsVec = [2, 10, 20, 50];

A = length(attackSampsVec);
S = length(sustainSampsVec);
R = length(releaseSampsVec);

Nfft = 2^nextpow2(N);
faxis = (fs/2) * linspace(0, 1, Nfft/2+1);

for i = 1:A    
    for j=1:S
        for k=1:R
            % excitation
            e = zeros(N,1); 
            
            attackSamps = attackSampsVec(i);
            sustainSamps = sustainSampsVec(j);
            releaseSamps = releaseSampsVec(k);
            
            e(1:attackSamps) = linspace(0, 1, attackSamps);
            e(attackSamps+1:attackSamps+1+sustainSamps) = 1;
            e(attackSamps+sustainSamps+1:attackSamps+sustainSamps+releaseSamps) = linspace(1, 0, releaseSamps);

            de = [diff(e); 0]; % take derivative for velocity?
            y = percSynthExcitation(de, yMSWav);
            y = y(1:N);
            
            y = y(attackSamps+sustainSamps+releaseSamps:end);
            
            if WRITETOFILE == 1
                audiowrite([outputDir 'ASR/' wavFileName(1:end-4) '_A' num2str(attackSamps) 'S' num2str(sustainSamps) 'R' num2str(releaseSamps) '.wav'], ...
                    scaleForSavingAudio(y), fs);
            end
            
            fprintf('attackSamps: %f, sustainSamps: %f, releaseSamps: %f\n', attackSamps, sustainSamps, releaseSamps);
            DE = fft(de, Nfft);
            DEPos = DE(1:Nfft/2+1);
            plot(faxis, abs(DEPos));

            
%         E = fft(e, Nfft);
%         EPos = E(1:Nfft/2+1);
%         plot(faxis, abs(EPos));
            soundsc(y,fs)
            keyboard

        end
    end
end


%% NOISE BURST

durNBVec = 0.01:0.05:0.2;
lowFreq = 20:200:600;
highFreq = [20000, 2000, 1000, 500];

D = length(durNBVec);
L = length(lowFreq);
H = length(highFreq);

 for j=1:L
        for k=1:H
               
                    for i=1:D

            
            % excitation
            e = zeros(N,1); 
            
            durNB = durNBVec(i);
            flow = lowFreq(j);
            fhigh = highFreq(k);
            
            % noise burst
            sampNB = ceil(durNB*fs);
            e = [(2*rand(1, sampNB))-1 zeros(1, N-sampNB)]';
            
            % filter
            [B, A] = butter(2, [flow/(fs/2) fhigh/(fs/2)], 'bandpass');
            e = filter(B, A, e);

            de = [diff(e); 0];
            y = percSynthExcitation(de, yMSWav);
            y = y(1:N);
            
            y = y(floor(sampNB/2):end);
            
            if WRITETOFILE == 1
                audiowrite([outputDir 'noiseBurst/' wavFileName(1:end-4) 'durNB' num2str(durNB) 'flow' num2str(flow) 'fhigh' num2str(fhigh) '.wav'], ...
                    scaleForSavingAudio(y), fs);
            end
            
            fprintf('durNB: %f, flow: %f, fhigh: %f\n', durNB, flow, fhigh);
            
            soundsc(y,fs)
            DE = fft(de, Nfft);
            DEPos = DE(1:Nfft/2+1);
            plot(faxis, abs(DEPos));
            keyboard
        end
    end
end


