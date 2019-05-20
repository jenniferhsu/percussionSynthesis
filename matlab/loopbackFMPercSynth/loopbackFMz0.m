function [z0, pitchGlide] = loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, varargin)
%%loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, varargin) generates 
% a loopback FM signal using the closed form loopback FM equation:
%   z_0(n) = (b0 + exp(j*w_0*n*T)) / (1 + b0 * exp(j*w_0*n*T))
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
%   b0: loopback FM closed-form timbre coefficient between -1 and 1
%       b0 is a single number if there is no time-varying timbre and
%       if there is timbre variation, b0 is a vector
%   dur: duration of sample in seconds
%   fs: sampling rate in Hz
%   varargin: optional structure that holds zc arguments for when 
%       pitchGlideType is 'exp', 'linB' or 'expB'
%       if pitchGlideType is 'exp', varargin holds values for T60
%       if pitchGlideType is 'linB', varargin holds values for
%           B, BEnd, and wc
%       if pitchGlideType is 'expB', varargin holds values for
%           g and wc
% outputs:
%   z0: the closed-form loopback FM signal
%   pitchGlide: the pitch glide 
%
% see tests/loopbackFMz0_Tests.m for examples
%
% author: Jennifer Hsu
% date: Spring 2019

N = fs*dur;                 % length of signal in samples
T = 1/fs;                   % sample period
n = linspace(0, N-1, N);    % sample vector

% if b0 is a single number, extend it to length N
if length(b0) == 1
    b0 = b0 * ones(1, N);
end

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
    a = varargin{1};
    T60 = a.T60;
    n60 = T60*fs;
    tau_d = -N*T/log(0.001/1);
    d = 1*exp(-n60*T/tau_d);
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
z0 = (b0 + exp(1j*Theta0)) ./ (1 + b0 .* exp(1j*Theta0));


end