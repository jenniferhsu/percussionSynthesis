function [yMS, yMSMat] = modalSynthSine(fVec, env, fs, pitchGlide)
%MODALSYNTHSINE.m performs modal synthesis with sinusoids instead of
%bandpass filters.
%
% inputs:
%   fVec: a vector of center frequencies in Hz (wc = 2*pi*fVec)
%       OR for a pitch glide, fVec is a matrix where size(fVec,1) = number
%       of center frequencies and every value in fVec is the desired pitch
%       glide value at a certain point in time.
%   env: amplitude envelope with length equal to the desired length for the
%       output signal.
%       ** If this is a matrix, it should be of dimensions
%       "length of fVec" by "desired output signal length" and each
%       envelope will be used for a different frequency in fVec.
%   fs: sampling rate in Hz
%   pitchGlide: set to 1 if we are making a pitch glide (don't multiply by
%       nT in the loop)
% outputs:
%   yMS: the output signal
%   yMSMat: the modal synthesis signal synthesized for each frequency in fVec
%       in matrix form. yFBFMMat(x,:) is the feedback FM signal created for
%       center frequency fVec(x)
%
% example:
% fs = 44100;
% dur = 1;
% fVec = 10^3 * [5.5000 7.6389 9.7778]; 
% env = g.^(linspace(0, N, N));
% [yMS, yMSMat] = modalSynthSine(fVec, env, fs);

% make sure parameters are all set
Nf = size(fVec, 1);
N = size(env, 2);
T = 1/fs;
nT = linspace(0, N*T, N);

% set up the envelope matrix
if size(env, 1)==1
    % repeat env as a matrix
    env = repmat(env, [Nf, 1]);
end

if size(fVec, 1) == 1 && size(fVec, 2) == 1
    % repeat fVec as a matrix
    fMat = repmat(fVec, [N, 1]);
elseif size(fVec, 1) == Nf && size(fVec, 2) == 1
    fMat = repmat(fVec, [N, 1]);
else
    fMat = fVec;
end

% synthesize modal synthesis signal (using sinusoids)
yMSMat = zeros(Nf, N); 

if pitchGlide == 0
    for i=1:Nf
        yMSMat(i,:) = yMSMat(i,:) + (exp(1j*2*pi*fMat(i,:).*nT) .* env(i,:));
    end
else
    % pitch glide
    for i=1:Nf
        yMSMat(i,:) = yMSMat(i,:) + (exp(1j*2*pi*fMat(i,:)) .* env(i,:));
    end
end

yMS = sum(yMSMat, 1);