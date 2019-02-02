function [y, yMat, TVAPFParamsUsed] = stretchedAPFAndTVAPFSynthesis(fVec, b0, env, TVAPFParams, r, fs, fVecEnd, pitchGlideMode)
%stretchedAPFAndTVAPFSynthesis generates a feedback FM signal using the
%stretched allpass filter and runs that through a 2nd order time-varying
%APF
%
% inputs:
%   fVec: a vector of sounding frequencies in Hz (w0 = 2*pi*fVec)
%       if pitchGlideMode is set to 'exp' or 'linear', this would be a
%       vector of starting sounding frequencies
%   b0: timbre control, I usually set this using B from FBFM rotation
%       equation
%   env: amplitude envelope with length equal to the desired length for the
%       output signal. 
%       ** If this is a matrix, it should be of dimensions
%       "length of fVec" by "desired output signal length" and each
%       envelope will be used for a different frequency in fVec.
%   TVAPFParams: struct that holds time-varying allpass filter parameters 
%       M, f_m, and f_b. 
%   r: if r==0, use fixed parameters for TVAPFParams. if r==1, use randomly
%       generated numbers for TVAPFParams
%   fs: sampling rate in Hz
%   fVecEnd: a vector of ending sounding frequencies in Hz with length
%       matching fVec. this can be an empty vector if pitchGlideMode is set
%       to 'none'
%   pitchGlideMode: for different pitch glides, set to
%       'none': no pitch glide
%       'exp': exponential pitch glide
%       'linear': linear pitch glide
%       'fbfm': rotation formulation type pitch glide with B=g^n. fVec is 
%           ignored and the starting frequency is 0 Hz. b0 is set according 
%           to B so that we also have a change in timbre over time
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
% fVecEnd = 4*fVec;
% B = 0.9;
% b0 = (sqrt(1-B^2) - 1)/B;
% env = g.^(linspace(0, N, N));
% TVAPFParams.M = fs/40;
% TVAPFParams.f_m = 100;
% TVAPFParams.f_b = fs/16;
% r = 1
% [y, yMat, TVAPFParamsUsed] = stretchedAPFAndTVAPFSynthesis(fVec, b0, env, TVAPFParams, r, fs, fVecEnd, 'fbfm');

% make sure parameters are all set
Nf = length(fVec);
N = size(env, 2);

dur = N/fs;
T = 1/fs;
nT = 0:T:((N/fs)-T);

% set up the envelope matrix
if size(env, 1)==1
    % repeat env as a matrix
    env = repmat(env, [Nf, 1]);
end

% TVAPFParameters
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

% synthesize time-varying APF and FBFM additive synthesis signal
yMat = zeros(Nf, N);

for i=1:Nf
    
    % synthesize sinusoid at correct frequency using FBFM stretched allpass
    % filter formulation
    f = fVec(i);
    
    % pitch glide modes
    if strcmp(pitchGlideMode, 'none')
        w0 = 2*pi*f*ones(1, N);
        b0Vec = b0 * ones(1, N);
        Theta = w0.*nT;
    elseif strcmp(pitchGlideMode, 'linear')
        m = (fVecEnd(i) - f)/(N/fs);
        b = f;
        %f0 = m*nT + b;
        %w0 = 2*pi*f0;
        
        b0Vec = b0 * ones(1, N);
        
        C = 0; 
        Theta = (m/2*nT.^2) + b*nT + C;
    elseif strcmp(pitchGlideMode, 'exp')
        tau = 1/(log(fVecEnd(i)/f));
        %f0 = f*exp(nT./tau);
        %w0 = 2*pi*f0;
        
        b0Vec = b0 * ones(1, N);
        
        C = 0;                      
        Theta = 2*pi*f*tau.*exp(nT./tau) + C;
    elseif strcmp(pitchGlideMode, 'fbfm')
        % fVec is ignored here and the starting freq is 0Hz
        % g can range from 0.9997 to 0.9999999, but it starts to break
        % down beyond that range
        f = fVecEnd(i);
        g = 0.9999;
        B = g.^(0:N-1);
        u = sqrt(1 - B.^2);
        
        b0Vec = (sqrt(1 - B.^2) - 1)./B;
        
        C = 0;                      % constant of integration
        Theta = (fVecEnd(i)*T/log(g)) .* (u - atanh(u) + C);
    end

    % stretched allpass filter synthesis
    x = (b0Vec + exp(1j*Theta)) ./ (1 + b0Vec.*exp(1j*Theta));
    if strcmp(pitchGlideMode, 'fbfm')
        x(1) = 0; % it computes as NAN 
    end
    
    % pass the oscillator through the time-varying APF
    %f_pi = fVec(i);
    f_pi = f;
    x = TVAPF2(x, f_pi, TVAPFParamsUsed.f_b(i), ...
                        TVAPFParamsUsed.M(i), ...
                        TVAPFParamsUsed.f_m(i), fs);
    
    % multiply with amplitude envelope
    yMat(i,:) = x .* env(i,:);
    
end
y = sum(yMat, 1);

end

