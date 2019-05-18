function [m, mMat] = loopbackFMMS(oscType, env, argStruct, fs)
% loopbackFMMS(oscType, env, argStruct, fs) creates a MS signal using
% sinusoidal, loopback FM zc, or looback FM z0 oscillators.  
%
% inputs:
%   oscType: oscillator type
%       's' - sinusoidal
%       'zc' - sample-by-sample loopback FM 
%       'z0' - closed-form loopback FM 
%   env: envelopes to use for each frequency component of size [Nf,N]
%       where Nf = number of modal frequencies and N = length of signal in
%       samples
%   argStruct: structure that holds arguments for each of the three
%       different types of oscillators. Only the struct for the oscillator
%       type that you are using must have values. See the struct
%       organization in the code below or look at examples in
%       tests/loopbackFMMS_Tests.m
%   fs: sample rate
%
% outputs:
%   m: the modal synthesis output signal
%   mMat: matrix holding output of each frequency component before
%       enveloping
%
% author: Jennifer Hsu
% date: Spring 2019

Nf = size(env,1);       % number of frequency components
N = size(env,2);        % length of output signal, samples
dur = N/fs;             % length of output signal, seconds

m = zeros(1, N);        % signal output
mMat = zeros(Nf, N);    % matrix of each frequency components output

if strcmp(oscType, 's')
    
    % traditional MS with sinusoidal oscillators
    f0Vec = argStruct.sinusoidArgs.f0Vec;
    f0EndVec = argStruct.sinusoidArgs.f0EndVec;
    pitchGlideTypeVec = argStruct.sinusoidArgs.pitchGlideTypeVec;
    zcArgsVec = argStruct.sinusoidArgs.zcArgsVec;
    
    for f=1:Nf
        f0 = f0Vec(f);
        f0End = f0EndVec(f);
        pitchGlideType = pitchGlideTypeVec{f};
        zcArgs = zcArgsVec(f);
        [mMat(f,:), ~] = sinusoid(f0, f0End, pitchGlideType, dur, fs, zcArgs);
        m = m + (mMat(f,:) .* env(f,:));
    end
    
elseif strcmp(oscType, 'zc')
    
    % sample-by-sample rotation loopback FM
    fcVec = argStruct.zcArgs.fcVec;
    BVec = argStruct.zcArgs.BVec;
    BEndVec = argStruct.zcArgs.BEndVec;
    gVec = argStruct.zcArgs.gVec;
    BGlideTypeVec = argStruct.zcArgs.BGlideType;
    
    for f=1:Nf
        fc = fcVec(f);
        B = BVec(f);
        BEnd = BEndVec(f);
        g = gVec(f);
        BGlideType = BGlideTypeVec(f);
        [mMat(f,:), ~] = loopbackFMzc(fc, B, BEnd, g, BGlideType, dur, fs);
        m = m + (mMat(f,:) .* env(f,:));
    end
    
elseif strcmp(oscType, 'z0')
    
    % closed-form loopback FM
    f0Vec = argStruct.z0Args.f0Vec;
    f0EndVec = argStruct.z0Args.f0EndVec;
    pitchGlideTypeVec = argStruct.z0Args.pitchGlideTypeVec;
    b0Mat = argStruct.z0Args.b0Mat;
    zcArgsVec = argStruct.z0Args.zcArgsVec;
    for f=1:Nf
        f0 = f0Vec(f);
        f0End = f0EndVec(f);
        pitchGlideType = pitchGlideTypeVec(f);
        b0 = b0Mat(f,:);
        zcArgs = zcArgsVec(f);
        [mMat(f,:), ~] = loopbackFMz0(f0, f0End, pitchGlideType, b0, dur, fs, zcArgs);
        m = m + (mMat(f,:) .* env(f,:));
    end
    
end

end