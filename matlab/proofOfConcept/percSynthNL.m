function y = percSynthNL(excitation, yMSWav, resIRwav)
%%Y = PERCSYNTHNL(EXCITATION, YMSWAV, RESIRWAV)
% This function uses modal synthesis percussion signals that
% have been modified with feedback FM and/or time-varying allpass filters
% and convolves them with excitation and body resonator signals.
% An excitation is convolved with the modal synthesis percussion signal.
% That aggregate signal is then convolved with a resonator impulse 
% response.  
% This is a test to see if there is a difference between this and
% percSynth() since my modal synthesis percussion signals are nonlinear.
%
% inputs:
%   excitation: the excitation function to use
%   yMSWav: the location of the modal synthesis percussion signal
%   resIRWav: the location of the body resonator impulse response
% outputs:
%   y: the synthesized percussion signal whose length will match the input
%       excitation
%
% see percSynthTestsNL.m for how to use this function to compare with
% results from percSynth()
% 

    N = length(excitation);

    % resonating body impulse response
    [resIR, ~] = audioread(resIRwav);
    resIR = resIR(:,1);

    % modal synthesis percussion sound
    [yMS, ~] = audioread(yMSWav);

    %% perform the convolutions

    % excitation convolved with modal synthesis result
    Nfft = 2^nextpow2(N + length(yMS) - 1);
    E = fft(excitation, Nfft);
    M = fft(yMS, Nfft);
    EM = E.*M;
    em = ifft(EM, Nfft);

    % aggregate excitation/modalsynthesis convolved % with resonating body
    Nfft = 2^nextpow2(length(em) + length(resIR) - 1);
    EM = fft(em, Nfft);
    R = fft(resIR, Nfft);
    Y = EM.*R;
    y = real(ifft(Y, Nfft));

end

