function [yFBFM, yFBFMMat] = feedbackFMSynthesis(fVec, BVec, env, fs)
%feedbackFMSynthesis generates a feedback FM signal using additive
% synthesis
%
% inputs:
%   fVec: a vector of center frequencies in Hz
%   BVec: for pitch glides, this is a vector of B coefficient values
%       if a pitch glide is not desired, this can be a single number
%   env: amplitude envelope with length equal to the desired length for the
%       output signal. 
%   fs: sampling rate in Hz
% outputs:
%   yFBFM: the output signal
%   yFBFMMat: the feedback FM signal synthesized for each frequency in fVec
%       in matrix form. yFBFMMat(x,:) is the feedback FM signal created for
%       center frequency fVec(x)
%
% example using the von Karman 3x3 steel plate modal frequencies:
% fs = 44100;
% dur = 1;
% g = 0.9999;
% N = fs*dur;
% fVec = 10^3 * [5.5000 7.6389 9.7778]; 
% BVec = 0.9; % or BVec = g.^(0:N-1)' for a pitch glide
% env = g.^(linspace(0, N, N));
% [yFBFM, yFBFMMat] = feedbackFMSynthesis(fVec, BVec, env, fs);

% make sure parameters are all set
N = length(env);
Nf = length(fVec);
T = 1/fs;

if length(BVec) == 1
    BVec = BVec * ones(1, N);
end

% synthesize feedback FM signal
yFBFMMat = zeros(Nf, N);
yFBFMMat(:,1) = 1;

for i=1:Nf
    f = fVec(i);
    wc = 2*pi*f;
    for n=2:N
        yFBFMMat(i,n) = exp(1j*wc*T*(1 + BVec(n-1) * real(yFBFMMat(i,n-1)))) * yFBFMMat(i,n-1);
    end
    yFBFMMat(i,:) = yFBFMMat(i,:).*env;
end
yFBFM = sum(yFBFMMat, 1);

end

