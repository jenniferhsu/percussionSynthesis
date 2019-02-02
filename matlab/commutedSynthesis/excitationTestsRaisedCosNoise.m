% percSynthTestsSMC2019.m
% this script shows how to use the percSynth function with different
% impulses

addpath('../proofOfConcept');

%% input parameters
fs = 44100;
N = 2*fs;

WRITETOFILE = 1; % set to 1 to save audio examples
outputDir = 'audioExamples/excitations/';
plotExcitations = 1;

% modal synthesis signal
%yMSWav = '../proofOfConcept/audioExamples/timeVaryingAPF/yTVAPFModal2.wav';
%yMSWav = '../proofOfConcept/audioExamples/FBFMStretchedAPF/ySAPFModal1.wav';
yMSWav = '../proofOfConcept/audioExamples/membrane/kettledrum-ySAPFMemb2.wav';

% input percussion signal for comparison
[yMS, ~] = audioread(yMSWav);

% get things ready for saving
if ~exist([outputDir 'raisedCosine/'])
    mkdir([outputDir 'raisedCosine/'])
end
if ~exist([outputDir 'noiseBurst/'])
    mkdir([outputDir 'noiseBurst/'])
end

tokens = strsplit(yMSWav, '/');
wavFileName = tokens{end};

%% Plot excitations
% this will stop running code after plotting the excitations

if plotExcitations==1
    
    % raised cosine
    winLength = 1024;
    e_rc = zeros(N,1); 
    n = winLength/2:winLength-1;
    w = 0.5 * (1 - cos((2*pi*n)/(winLength-1)));

    e_rc(1:winLength/2) = ones(winLength/2, 1);
    e_rc(winLength/2+1:winLength) = w;
    
    % filtered noise burst
    e_nb = zeros(N,1); 

    durNB = 0.03;
    flow = 200;
    fhigh = 2000;

    % noise burst
    sampNB = ceil(durNB*fs);
    e_nb = [(2*rand(1, sampNB))-1 zeros(1, N-sampNB)]';

    % filter
    [B, A] = butter(2, [flow/(fs/2) fhigh/(fs/2)], 'bandpass');
    e_nb = filter(B, A, e_nb);
    
    figure
    
    subplot(211)
    plot(e_rc, 'linewidth', 2);
    xlabel('Time (samples)');
    ylabel('Amplitude (linear)');
    title('Raised cosine excitation, window length L=1024 samples');
    xlim([0 2048]);
    ylim([0 1.1]);
    grid on
    set(gca, 'FontSize', 15);
    
    subplot(212)
    plot(e_nb, 'linewidth', 2);
    xlabel('Time (samples)');
    ylabel('Amplitude (linear)');
    title('Filtered noise burst excitation, \tau _d = 0.01 sec, fc_{low} = 200, fc_{high} = 2000');
    xlim([0 2048]);
    ylim([-0.8 0.8]);
    grid on
    set(gca, 'FontSize', 15);
    
    saveas(gcf, 'figures/excitations', 'epsc');
    
    return
end


%% RAISED COSINE
% making winLength shorter retains more of the higher frequencies
% longer winLength removes the higher frequencies and gives a bassier sound
% the curve doesn't make a difference, so let's change it to the end of a
% cosine

winLengthVec = [8, 16, 64, 256, 1024];
Nwl = length(winLengthVec);

Nfft = 2^nextpow2(N);
faxis = (fs/2) * linspace(0, 1, Nfft/2+1);

for i = 1:Nwl   

    % excitation
    e = zeros(N,1); 

    winLength = winLengthVec(i);

    % raised cosine/Hann window
    n = winLength/2:winLength-1;
    w = 0.5 * (1 - cos((2*pi*n)/(winLength-1)));

    e(1:winLength/2) = ones(winLength/2, 1);
    e(winLength/2+1:winLength) = w;
    
    % take derivative for velocity
    de = [diff(e); 0]; 
    y = percSynthExcitation(de, yMSWav);
    y = y(1:N);

    %y = y(winLength:end);

    if WRITETOFILE == 1
        audiowrite([outputDir 'raisedCosine/' wavFileName(1:end-4) '_winLength_' num2str(winLength) '_P1x_' num2str(P1x) '_P1y_' num2str(P1y) '.wav'], ...
            scaleForSavingAudio(y), fs);
    end

    fprintf('winLength: %f, P1x: %f, P1y: %f\n', winLength, P1x, P1y);
%     DE = fft(de, Nfft);
%     DEPos = DE(1:Nfft/2+1);
%     plot(faxis, abs(DEPos));
% 
% 
%     E = fft(e, Nfft);
%     EPos = E(1:Nfft/2+1);
%     plot(faxis, abs(EPos));
%     soundsc(y,fs)
%     keyboard

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
             
%             soundsc(y,fs)
%             DE = fft(de, Nfft);
%             DEPos = DE(1:Nfft/2+1);
%             plot(faxis, abs(DEPos));
%             keyboard
        end
    end
end


