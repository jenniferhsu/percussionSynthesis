function [s, pitchGlide] = sinusoid(f0, f0End, pitchGlideType, dur, fs, varargin)
%sinusoid(f0, b0, dur, fs) generates a sinusoidal signal for traditional
%   modal synthesis with the ability to create pitch glides
%
% inputs:
%   f0: starting sounding frequency in Hz (w0 = 2*pi*f0)
%       f0 is a single number whether or not there is a pitch glide
%   f0End: ending sounding frequency in Hz (w0End = 2*pi*f0End)
%       f0End is a single number whether or not there is a pitch glide
%   pitchGlideType: 'none', 'lin', 'exp', 'sqrt', 'linB', 'expB' 
%       'lin', 'exp', and 'sqrt' specify functions that define the pitch 
%           trajectory from f0 to f0End
%       'linB' and 'expB' mean that feedback coefficient B from the 
%           sample-by-sample rotation equation (z_c(n)) moves linearly 
%           and exponentially from B to f0 to B = f0End   
%   dur: duration of sample in seconds
%   fs: sampling rate in Hz
%   varargin: optional structure that holds zc arguments for when 
%       pitchGlideType is 'linB' or 'expB'
%       if pitchGlideType is 'linB', varargin holds values for
%           B, BEnd, and wc
%       if pitchGlideType is 'expB', varargin holds values for
%           g and wc
% outputs:
%   s: the resulting sinusoid
%   pitchGlide: the pitch glide 
%
% see tests/sinusoid_Tests.m for examples
%
% author: Jennifer Hsu
% date: Spring 2019

N = fs*dur;                 % length of signal in samples
T = 1/fs;                   % sample period
n = linspace(0, N-1, N);    % sample vector

% figure out pitch glide function
Theta0 = zeros(1, N);
if strcmp(pitchGlideType, 'none')
    Theta0 = 2*pi*f0*n*T;
    pitchGlide = [];
elseif strcmp(pitchGlideType, 'lin')
    k = (f0End - f0)/dur;
    l = f0;
    pitchGlide = k*n*T + l;
    Theta0 = 2*pi*((k*(n*T).^2)/2 + l*n*T);
elseif strcmp(pitchGlideType, 'exp')
    tau_d = -N*T/log(0.001/1);
    d = 1*exp(-n*T/tau_d);
    pitchGlide = (f0 - f0End)*d + f0End;
    Theta0 = -2*pi*(f0 - f0End)*tau_d*1*exp(-n*T/tau_d) + 2*pi*f0End*n*T;
elseif strcmp(pitchGlideType, 'sqrt')
    a = (f0End - f0)/sqrt(N*T);
    c = f0;
    pitchGlide = a*sqrt(n*T) + c;
    Theta0 = 2*pi*((2/3)*a*(n*T).^(3/2) + c*n*T);
elseif strcmp(pitchGlideType, 'linB')
    zcArgs = varargin{1};
    B = zcArgs.B;
    BEnd = zcArgs.BEnd;
    wc = zcArgs.wc;
    k = (BEnd - B)/N;
    l = B;
    BVec = k*n + l;
    w0 = wc * sqrt(1 - BVec.^2);
    pitchGlide = w0/(2*pi);
    Theta0 = (wc*T/(2*k)) .* ...
        ((k*n + l) .* sqrt(1 - (k*n + l).^2) + asin(k*n + l));  
elseif strcmp(pitchGlideType, 'expB')
    zcArgs = varargin{1};
    g = zcArgs.g;
    wc = zcArgs.wc;
    BVec = g.^n;
    w0 = wc * sqrt(1 - BVec.^2);
    pitchGlide = w0/(2*pi);
    Theta0 = wc*T/log(g) .* (sqrt(1 - g.^(2*n)) - atanh(1 - g.^(2*n)));
end
    

% the closed-form loopback FM calculation
s = exp(j*Theta0);


end