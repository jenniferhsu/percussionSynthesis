function [ySAPF, ySAPFMat] = stretchedAPFSynthesis(fVec, b0, env, fs, fVecEnd, pitchGlideMode)
%stretchedAPFSynthesis generates a feedback FM signal using additive
% synthesis with the stretched APF formulation of FBFM
%
% inputs:
%   fVec: a vector of sounding frequencies in Hz (w0 = 2*pi*fVec)
%       if pitchGlideMode is set to 'exp' or 'linear', this would be a
%       vector of starting sounding frequencies
%   b0: timbre control, I usually set this using B from FBFM rotation
%       equation
%   env: amplitude envelope with length equal to the desired length for the
%       output signal. 
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
%   ySAPF: the output signal
%   ySAPFMat: the feedback FM signal synthesized for each frequency in fVec
%       in matrix form. yFBFMMat(x,:) is the feedback FM signal created for
%       center frequency fVec(x)
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
% [ySAPF, ySAPFMat] = stretchedAPFSynthesis(fVec, b0, env, fs, fVecEnd, 'linear');

% set this to 1 if you want to make plots for SMC 2019 paper
SMC2019PLOT = 0;

% make sure parameters are all set
N = length(env);
Nf = length(fVec);
T = 1/fs;
nT = 0:T:((N/fs)-T);

% synthesize feedback FM signal using the strecthed allpass filter
ySAPFMat = zeros(Nf, N);

% for SMC 2019 plotting
osc = zeros(Nf, N);

for i=1:Nf
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
        g = 0.9999;
        B = g.^(0:N-1)';
        u = sqrt(1 - B.^2);
        
        b0Vec = (sqrt(1 - B.^2) - 1)./B;
        
        C = 0;                      % constant of integration
        Theta = (fVecEnd(i)*T/log(g)) .* (u - atanh(u) + C);
    end
    
    % stretched allpass filter
    ySAPFMat(i,:) = (b0Vec + exp(1j*Theta)) ./ (1 + b0Vec.*exp(1j*Theta));
    osc(i,:) = ySAPFMat(i,:); % for SMC 2019 plotting
     
%     relationships that might be useful
%     wc = w0 / sqrt(1 - ((-2*b0)/(b0^2 + 1))^2);
%     wc = 2*pi*f;
%     w0 = wc*sqrt(1 - BVec(1)^2);
%     b0 = ((w0./wc) - 1) ./ sqrt(1 - w0.^2./wc^2);
   
%     figure
%     spectrogram(real(ySAPFMat(i,:)), hann(256), 128, 1024, fs, 'yaxis');

%     fourier transform check of this
%     Nfft = 2^nextpow2(N);
%     faxis = (fs/2)*linspace(0, 1, Nfft/2+1);
%     
%     YFBFM = fft(yFBFMMat(i,:), Nfft);
%     YFBFMPos = YFBFM(1:Nfft/2+1);
%     
%     YSAPF = fft(ySAPFMat(i,:), Nfft);
%     YSAPFPos = YSAPF(1:Nfft/2+1);
%     
%     figure
%     plot(faxis, abs(YFBFMPos))
%     hold on
%     plot(faxis, abs(YSAPFPos), 'r--')
%     plot([w0/(2*pi) w0/(2*pi)], [0 max(abs(YFBFMPos))], 'k--');
%     keyboard

    ySAPFMat(i,:) = ySAPFMat(i,:).*env;
end
ySAPF = sum(ySAPFMat, 1);

if SMC2019PLOT==1
    
    % LOOPBACK FM WITH ADDITIVE SYNTHESIS
    figure
    subplot(311); plot(real(osc(1, 1:100)), 'linewidth', 2)
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    
    subplot(312); plot(real(osc(2, 1:100)), 'linewidth', 2)
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    
    subplot(313); plot(real(osc(3, 1:100)), 'linewidth', 2)
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    print(['~/Documents/ucsd/winter2019/smc/background_planning/figure_making/' 'loopbackFMOsc'], '-deps', '-r0')

    % ENVELOPE
    figure
    subplot(311); plot(real(env), 'linewidth', 2)
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    
    subplot(312); plot(real(env), 'linewidth', 2)
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    
    subplot(313); plot(real(env), 'linewidth', 2)
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 6];
    print(['~/Documents/ucsd/winter2019/smc/background_planning/figure_making/' 'loopbackFMEnv'], '-deps', '-r0')
    
    
    % ADDITIVE SYNTHESIS RESULT
    figure
    plot(nT, ySAPF, 'linewidth', 2);
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 3];
    print(['~/Documents/ucsd/winter2019/smc/background_planning/figure_making/' 'additiveSynthesisResult'], '-deps', '-r0')
keyboard
end

end

