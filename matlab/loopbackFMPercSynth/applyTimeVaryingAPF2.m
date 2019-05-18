function [y, yMat] = applyTimeVaryingAPF2(mMat, env, fs, TVAPFParams)
%%applyTimeVaryingAPF2(mMat, env, fs, TVAPFParams) applies time-varying
% 2nd order allpass filters to a matrix of loopback FM (or traditional MS)
% signals.  These matrices are the ones returned by loopbackFMMS.
%
% inputs:
%   mMat: matrix of loopback FM or traditional MS signals as returned by
%       the loopbackFMMS function. the size is [Nf, N] where
%       Nf = number of modal components
%       N = lenght of signal in samples
%   env: amplitude envelope matrix of size [Nf, N]
%   fs: sampling rate (samples/sec)
%   TVAPFParams: struct that holds time-varying allpass filter parameters 
%       fpi, M, f_m, and f_b as vectors of length Nf
%
% output:
%   y: the time-varying, 2nd-order allpassed loopback FM signal
%   yMat: a matrix holding the individual time-varying, 2nd-order 
%       allpassed loopback FM oscillators before enveloping
%
% see tests/applyTimeVaryingAPF2_Tests.m for examples

Nf = size(mMat, 1);
N = size(mMat, 2);

fpiVec = TVAPFParams.fpiVec;
fbVec = TVAPFParams.fbVec;
MVec = TVAPFParams.MVec;
fmVec = TVAPFParams.fmVec;

yMat = zeros(Nf, N);
y = zeros(1, N);

for f=1:Nf
    fpi = fpiVec(f);
    fb = fbVec(f);
    M = MVec(f);
    fm = fmVec(f);
    yMat(f,:) = timeVaryingAPF2(real(mMat(f,:)), fpi, fb, M, fm, fs);
    y = y + (yMat(f,:) .* env(f,:));
end




end