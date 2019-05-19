function [zc, pitchGlide] = loopbackFMzc(fc, B, BEnd, g, BGlideType, dur, fs)
%loopbackFMzc(fc, B, BEnd, g, BGlideType, dur, fs) generates a loopback FM 
% signal using the sample-by-sample rotation equation:
%   z_c(n) = exp(j*w_c*T*(1 + B * real(z_c(n-1)))) z_c(n-1)
%
% inputs:
%   fc: center frequency in Hz (wc = 2*pi*fc)
%   B: loopback FM coefficient between -1 and 1
%       starting B value if B is time-varying
%   BEnd: loopback FM coefficient between -1 and 1
%       ending B value if B is time-varying
%   g: exponetial function bass if BGlideType is expB
%   BGlideType: 'none', 'useB', 'linB', or 'expB' 
%       if 'none': fc and B are used in the loopback FM equation
%       if 'useB': the pitch glide is stored in B and nothing needs
%           to be calculated
%       if 'linB': fc, B, and BEnd are used for the time-variation
%       if 'expB': fc and g used for the time-variation
%   dur: duration of sample in seconds
%   fs: sampling rate in Hz
% outputs:
%   zc: the loopback FM signal
%   pitchGlide: the pitch glide 
%
% see tests/loopbackFMzc_Tests.m for examples
%
% author: Jennifer Hsu
% date: Spring 2019

N = fs*dur;                 % length of signal in samples
T = 1/fs;                   % sample period
wc = 2*pi*fc;               % carrier frequency in rad/sec
nn = linspace(0, N-1, N);    % sample vector

zc = zeros(1, N);
zc(1) = 1;

if strcmp(BGlideType, 'none')
    % static pitch and timbre
    for n=2:N
        zc(n) = exp(1j*wc*T*(1 + B * real(zc(n-1)))) * zc(n-1);
    end
    pitchGlide = [];
else
    % time-varying pitch and timbre
    if strcmp(BGlideType, 'useB')
        BVec = B;
        pitchGlide = fc*sqrt(1 - BVec.^2);
    elseif strcmp(BGlideType, 'linB')  
        k = (BEnd - B)/N;
        l = B;
        BVec = k*nn + l;
        w0 = wc * sqrt(1 - BVec.^2);
        pitchGlide = w0/(2*pi);
     elseif strcmp(BGlideType, 'expB')
        BVec = g.^nn;
        w0 = wc * sqrt(1 - BVec.^2);
        pitchGlide = w0/(2*pi);
    end
    for n=2:N
        zc(n) = exp(1j*wc*T*(1 + BVec(n-1) * real(zc(n-1)))) * zc(n-1);
    end
end
    

end