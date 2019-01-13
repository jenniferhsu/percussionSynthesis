function y = percSynth(excitation, yMSWav, resIRwav)
%%Y = PERCSYNTH(EXCITATION, YMSWAV, RESIRWAV)
% This function combines the ideas in commuted synthesis
% and convolutional synthesis with modal synthesis percussion signals that
% have been modified with feedback FM and/or time-varying allpass filters.
% An excitation is convolved with a resonator impulse response.  Then that 
% aggregate excitation is convolved with the modal synthesis percussion
% signal.
%
% inputs:
%   excitation: the excitation function to use
%   yMSWav: the location of the modal synthesis percussion signal
%   resIRWav: the location of the body resonator impulse response
% outputs:
%   y: the synthesized percussion signal whose length will match the input
%       excitation
%
% see percSynthTests.m for how to use this function
% 

    N = length(excitation);

    % resonating body impulse response
    [resIR, ~] = audioread(resIRwav);
    resIR = resIR(:,1);

    % modal synthesis percussion sound
    [yMS, ~] = audioread(yMSWav);

    %% perform the convolutions

    % excitation convolved with resonating body to form aggregate excitation
    Nfft = 2^nextpow2(N + length(resIR) - 1);
    E = fft(excitation, Nfft);
    R = fft(resIR, Nfft);
    ER = E.*R;
    er = ifft(ER, Nfft);

    % aggregate excitation convolved with mesh impulse response
    Nfft = 2^nextpow2(length(er) + length(yMS) - 1);
    AE = fft(er, Nfft);
    M = fft(yMS, Nfft);
    Y = AE.*M;
    y = real(ifft(Y, Nfft));

end

