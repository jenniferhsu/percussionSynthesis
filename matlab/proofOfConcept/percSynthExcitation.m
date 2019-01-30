function y = percSynthExcitation(excitation, yMSWav)
%%Y = PERCSYNTHEXCITATION(EXCITATION, YMSWAV, RESIRWAV)
% This function uses modal synthesis percussion signals that
% have been modified with feedback FM and/or time-varying allpass filters
% and convolves them with an excitation signal.
%
% inputs:
%   excitation: the excitation function to use
%   yMSWav: the location of the modal synthesis percussion signal
% outputs:
%   y: the synthesized percussion signal whose length will match the input
%       excitation
%
% see percSynthTestsExcitation.m for how to use this function
% 

    N = length(excitation);

    % modal synthesis percussion sound
    [yMS, ~] = audioread(yMSWav);

    %% perform the convolutions

    % excitation convolved with modal synthesis result
    Nfft = 2^nextpow2(N + length(yMS) - 1);
    E = fft(excitation, Nfft);
    M = fft(yMS, Nfft);
    Y = E.*M;
    y = ifft(Y, Nfft);

end

