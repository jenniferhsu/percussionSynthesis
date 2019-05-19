function y = applyCommutedSynthesis(yMS, resIRWav, excitationType, excitationParams, fs)
%%applyCommutedSynthesis(yMS, resIRWav, excitationType, excitationParams,
% fs) performs commuted synthesis.
% inputs:
%   yMS: modal synthesis percussion sound
%   resIRWav: path to resonating IR audio file
%   excitationType: 'rc' for raised cosine and 'nb' for filtered 
%       noise burst
%   excitationParams: struct holding parameters for the excitation
%       if excitationType is 'rc':
%           excitationParams.winLength - window length
%       if excitationType is 'nb':
%           excitationParams.durNB - noise burst duration
%           excitationParams.lowFreq - low frequency cutoff for bp filter
%           excitationParams.highFreq - high frequency cutoff for bp filter
%   fs: sample rate (samples/sec)
% outputs:
%   y: the commuted synthesis output signal
%
%
% author: Jennifer Hsu
% date: Spring 2019

N = length(yMS);
excitation = zeros(N, 1);

%% create the excitation
if strcmp(excitationType, 'rc')
    % raised cosine/Hann window
    winLength = excitationParams.winLength;

    % calculate window
    n = [0:winLength-1]';
    excitation(1:winLength) = 0.5 * (1 - cos((2*pi*n)/(winLength-1)));

    % take first difference for "velocity"
    dexcitation = [diff(excitation); 0]; 

elseif strcmp(excitationType, 'nb')
    % noise burst excitation
    durNB = excitationParams.durNB;
    lowFreq = excitationParams.lowFreq;
    highFreq = excitationParams.highFreq;

    % calculate noise burst
    sampNB = ceil(durNB*fs);
    excitation(1:sampNB) = 2*rand(1, sampNB) - 1;%[(2*rand(1, sampNB))-1 zeros(1, N-sampNB)]';

    [B, A] = butter(5, [lowFreq/(fs/2) highFreq/(fs/2)], 'bandpass');
    %freqz(B,A)
    excitation = filter(B, A, excitation);

    % take derivative for velocity
    dexcitation = [diff(excitation); 0]; 
    
else
    disp('please check that your excitation type is "rc" or "nb"');
end


%% perform the commuted synthesis
N = length(excitation);

% resonating body impulse response
[resIR, ~] = audioread(resIRWav);
resIR = resIR(:,1);


%% perform the convolutions

% excitation convolved with resonating body to form aggregate excitation
Nfft = 2^nextpow2(N + length(resIR) - 1);
E = fft(dexcitation, Nfft);
R = fft(resIR, Nfft);
ER = E.*R;
er = ifft(ER, Nfft);

% aggregate excitation convolved with mesh impulse response
Nfft = 2^nextpow2(length(er) + length(yMS) - 1);
AE = fft(er, Nfft);
M = fft(yMS, Nfft);

if size(AE, 1) == 1
    if size(M,1) ~= 1
        M = M';
    end
elseif size(AE, 2) == 1
    if size(M,2) ~= 1
        M = M';
    end
end

Y = AE.*M;
y = real(ifft(Y, Nfft));

y = y(1:N);
