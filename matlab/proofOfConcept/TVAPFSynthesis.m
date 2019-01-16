function [yTVAPF, yTVAPFMat, TVAPFParamsUsed] = TVAPFSynthesis(fVec, env, TVAPFParams, r, fs)
%TVAPFSynthesis generates a feedback FM signal using additive
% synthesis
%
% inputs:
%   fVec: a vector of center frequencies in Hz
%   env: amplitude envelope with length equal to the desired length for the
%       output signal. 
%   TVAPFParams: struct that holds time-varying allpass filter parameters 
%       M, f_m, and f_b. 
%   r: if r==0, use fixed parameters for TVAPFParams. if r==1, use randomly
%       generated numbers for TVAPFParams
%   fs: sampling rate in Hz
% outputs:
%   yTVAPF: the output signal
%   yTVAPFMat: the feedback FM signal synthesized for each frequency in fVec
%       in matrix form. yFBFMMat(x,:) is the feedback FM signal created for
%       center frequency fVec(x)
%   TVAPFParamsUsed: time-varying allpass filter parameters that were
%       actually used for synthesis. this may differ from the input
%       TVAPFParams if r > 0
%
% example using the von Karman 3x3 steel plate modal frequencies:
% fs = 44100;
% dur = 1;
% g = 0.9999;
% N = fs*dur;
% fVec = 10^3 * [5.5000 7.6389 9.7778]; 
% env = g.^(linspace(0, N, N));
% TVAPFParams.M = fs/40;
% TVAPFParams.f_m = 100;
% TVAPFParams.f_b = fs/16;
% r = 1
% [yTVAPF, yTVAPFMat] = TVAPFSynthesis(fVec, env, fs);

% make sure parameters are all set
N = length(env);
dur = N/fs;
Nf = length(fVec);
T = 1/fs;
nVec = 0:T:(dur-T);

M = TVAPFParams.M;
f_m = TVAPFParams.f_m;
f_b = TVAPFParams.f_b;

if r == 0
    % f_b = 0 creates a clean sound
    % f_b = 1 creates a "WHIRRING", "WHIZZING" kind of sound 
    % as we go from f_b to 20000, the sound gradually becomes cleaner and purer
    % from f_b = 20000 up to fb = 198452, the sound stays clean
    % f_b > 198452 creates an unstable signal
    TVAPFParamsUsed.f_b = f_b * ones(Nf,1);

    % M can range from 0 to 2000000
    % the smaller M is, the less crazy it sounds
    % the larger M is, the noisier it sounds, "METALLIC WHITE NOISE"?
    % wow, this makes it sound really interesting. we might need to
    % make this it's own parameter
    TVAPFParamsUsed.M = M * ones(Nf,1); %2000000 * ones(Nf,1);

    % increasing f_m lowers the sounding frequency, "PITCH?"
    % f_m = 0 seems like the highest pitch
    % as f_m increases, the pitch lowers
    % when we get around 19000000, it starts to sound like two pitches
    % somewhere between 10^7 and 10^8, it jumps back from
    % lower to higher (some aliasing maybe?)
    TVAPFParamsUsed.f_m = f_m * ones(Nf,1); %19000000 * ones(Nf,1);
else
    max_M = 2000000;
    min_M = 0;
    TVAPFParamsUsed.M = max(0, (max_M-min_M).*rand(Nf,1) + min_M);

    max_f_m = 19000000;
    min_f_m = 0;
    TVAPFParamsUsed.f_m = max(0, (max_f_m-min_f_m).*rand(Nf,1) + min_f_m);

    max_f_b = 20000;
    min_f_b = 1;
    TVAPFParamsUsed.f_b = max(0, (max_f_b-min_f_b).*rand(Nf,1) + min_f_b);
end

% synthesize time-varying APF additive synthesis signal
yTVAPFMat = zeros(Nf, N);
yTVAPFMat(:,1) = 1;

for i=1:Nf
    
    % synthesize sinusoid at correct frequency
    f = fVec(i);
    x = exp(1j*2*pi*f*nVec);
    
    % pass the oscillator through the time-varying APF
    f_pi = fVec(i);
    x = TVAPF2(x, f_pi, TVAPFParamsUsed.f_b(i), ...
                        TVAPFParamsUsed.M(i), ...
                        TVAPFParamsUsed.f_m(i), fs);
    
    % multiply with amplitude envelope
    yTVAPFMat(i,:) = x .* env;
    
end
yTVAPF = sum(yTVAPFMat, 1);

end

