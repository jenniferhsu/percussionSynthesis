fs = 44100;
dur = 2;

addpath(genpath('../loopbackFMPercSynth/'));
addpath(genpath('../helperFunctions/'));
savePlots = 1;
saveAudio = 1;

figDir = 'figures/';
audioDir = 'audioExamples/';

% === derived parameters === 
% feedback FM
B = 0.3;    % feedback coefficient
g = 0.9999; % pitch glide coefficient

N = fs*dur;
nVec = 0:T:(dur-T);
T = 1/fs;
BVec = g.^(0:N-1);             % feedback FM pitch glide coefficient

% === steel pan modal frequencies === 
% from 
% https://courses.physics.illinois.edu/phys406/sp2017/Student_Projects/Spring16/Lienne_Pyzik_Physics_406_Final_Report_Sp16.pdf

fVec = [233, 468, 706, 1049, 1176, 2113];

% marimba
% f1 = 440;
% Nf = 7;
% fVec = zeros(1, Nf);
% for i=1:Nf
%     if i==1
%         fVec(i) = 3.011^2;
%     else
%         fVec(i) = (2*i+1)^2;
%     end
% end
% fVec = fVec/fVec(1);
% fVec = f1 * fVec;

Nf = length(fVec);

% center frequencies for feedback FM if sounding frequencies = membrane
% modal freqeuncies
fcVec = fVec./(sqrt(1-B^2));

% === set up decay envelopes ===
% want the lowest modal frequencies to have a larger starting amplitude 
% than the highest modal frequencies
aStart = zeros(1, Nf);
aStart(1) = 1;              % highest amplitude for lowest modal freq
aStart(Nf) = 0.01;          % lowest amplitude for highest modal freq
tau = -(Nf-1) / log(aStart(Nf));
a0 = 1/exp(-1/tau);
aStart(2:Nf-1) = a0*exp(-(2:Nf-1)/tau);     % exponentially decreasing 
                                            % starting amplitudes

e = g.^(linspace(0, N, N));   % exponential decay

env = zeros(Nf, N);
for i=1:Nf
    env(i,:) = aStart(i) * e;
end
wg = env(1,:)';

% === steel pan synthesis - loopback FM parameters ===

% traditional MS parameters
argStruct.sinusoidArgs.f0Vec = fVec;
argStruct.sinusoidArgs.f0EndVec = fVec;
argStruct.sinusoidArgs.pitchGlideTypeVec = cell(1, Nf);
for f=1:Nf
    argStruct.sinusoidArgs.pitchGlideTypeVec{f} = 'none';
end
argStruct.sinusoidArgs.zcArgsVec = zeros(1, Nf);

% loopback FM z0 parameters
argStruct.z0Args = struct();
argStruct.z0Args.f0Vec = fVec*1.05;
argStruct.z0Args.f0EndVec = fVec;
argStruct.z0Args.pitchGlideTypeVec = cell(1, Nf);
for f=1:Nf
    argStruct.z0Args.pitchGlideTypeVec{f} = 'exp';
    argStruct.z0Args.b0Mat(f,:) = linspace(0.3*((Nf-(f-1))/Nf), 0.001, N);
    argStruct.z0Args.zcArgsVec(f) = struct();
end
for f=1:Nf
    argStruct.z0Args.zcArgsVec(f).T60 = 0.9;
end

% === steel pan synthesis - traditional MS and loopback FM MS ===

% steel pan - traditional MS
[steelpan_s, steelpan_s_Mat] = loopbackFMMS('s', env, argStruct, fs);

% steel pan - loopback FM z0
[steelpan_z0, steelpan_z0_Mat] = loopbackFMMS('z0', env, argStruct, fs);


% === time-varying APF ===

% parameters
TVAPFParams.fpiVec = argStruct.sinusoidArgs.f0Vec;
TVAPFParams.MVec = [1000 2000 1500 3000 2300 6000 2500]/100;
TVAPFParams.fmVec = [78, 31, 42, 83, 100, 400, 300];
TVAPFParams.fbVec = [100 200 300 400 500 600 700]/10;

% time-varying APF - traditional MS
[ySteelpan_s, ySteelpan_s_Mat] = applyTimeVaryingAPF2(steelpan_s_Mat, env, fs, TVAPFParams);

% time-varying APF - loopback FM z0
[ySteelpan_z0, ySteelpan_z0_Mat] = applyTimeVaryingAPF2(steelpan_z0_Mat, env, fs, TVAPFParams);


% === steelpan synthesis - commuted synthesis ===

% commuted synthesis parameters
% RAISED COSINE
resIRWav = '../loopbackFMPercSynth/resonatorIRs/CarpenterCenter.wav';
excitationType = 'rc';
excitationParams = struct();
excitationParams.winLength = 8;

% steelpan - traditional MS commuted synthesis
steelpan_s_CS = applyCommutedSynthesis(real(steelpan_s), resIRWav, excitationType, excitationParams, fs);
%marimba_s_CS = marimba_s_CS .* wg;

% steelpan - loopback FM z0 commuted synthesis
steelpan_z0_CS_rc = applyCommutedSynthesis(real(steelpan_z0), resIRWav, excitationType, excitationParams, fs);

% NOISE BURST
resIRWav = '../loopbackFMPercSynth/resonatorIRs/CarpenterCenter.wav';
excitationType = 'nb';
excitationParams = struct();
excitationParams.durNB = .01;
excitationParams.lowFreq = 200;
excitationParams.highFreq = 8000;

steelpan_z0_CS_nb = applyCommutedSynthesis(real(steelpan_z0), resIRWav, excitationType, excitationParams, fs);


% NOISE BURST on TV AP
steelpan_z0_CS_nb_TVAP = applyCommutedSynthesis(real(ySteelpan_z0), resIRWav, excitationType, excitationParams, fs);


% === steelpan synthesis - commuted synthesis - processing a file ===

resIRWav = '../proofOfConcept/resonatorIRs/charm_08.wav';
excitationType = 'nb';
excitationParams = struct();
excitationParams.durNB = .01;
excitationParams.lowFreq = 200;
excitationParams.highFreq = 8000;

steelpan_z0_CS_nb_charm = applyCommutedSynthesis(real(steelpan_z0), resIRWav, excitationType, excitationParams, fs);


%% write audio files
if saveAudio
    audiowrite([audioDir 'steelpan_s' '.wav'], ...
                scaleForSavingAudio(real(steelpan_s)), fs);
    audiowrite([audioDir 'steelpan_z0' '.wav'], ...
                scaleForSavingAudio(real(steelpan_z0)), fs);
    audiowrite([audioDir 'steelpan_z0_TVAP' '.wav'], ...
                scaleForSavingAudio(real(ySteelpan_z0)), fs);
    audiowrite([audioDir 'steelpan_z0_commutedSynthesis_rc_CarpenterCenter' '.wav'], ...
                scaleForSavingAudio(real(steelpan_z0_CS_rc)), fs);
    audiowrite([audioDir 'steelpan_z0_commutedSynthesis_nb_CarpenterCenter' '.wav'], ...
                scaleForSavingAudio(real(steelpan_z0_CS_nb)), fs);
    audiowrite([audioDir 'steelpan_z0_commutedSynthesis_nb_CarpenterCenter_TVAP' '.wav'], ...
                scaleForSavingAudio(real(steelpan_z0_CS_nb_TVAP)), fs);
    audiowrite([audioDir 'steelpan_z0_commutedSynthesis_nb_charm' '.wav'], ...
                scaleForSavingAudio(real(steelpan_z0_CS_nb_charm)), fs);
            
    % save invidiual loopback FM envelopes oscillators
    for i=1:Nf
        audiowrite([audioDir 'steelpan_z0_' num2str(i) '.wav'], ...
                scaleForSavingAudio(env(1,:) .* real(steelpan_z0_Mat(i,:))), fs);
    end
end



%% plots

if savePlots
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% MODAL SYNTHESIS %%% 
    
    % sinusoids
    M = 175; % end index
    figure

    for i=1:Nf
        subplot(Nf, 1, i)
        plot(nVec(1:M), real(steelpan_s_Mat(i,1:M)), 'linewidth', 2);
        xlim([nVec(1) nVec(M)]);
        set(gca,'linewidth', 3)
        set(gca,'XTick',[], 'YTick', [])
    end

    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    print([figDir 'modalSynthesisDiagramSinusoids'], '-depsc', '-r0')
    
    % MS output
    figure
    plot(nVec, real(steelpan_s/max(abs(real(steelpan_s)))), 'linewidth', 2);
    xlim([nVec(1) nVec(end)]);
    ylim([-1 1]);
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([figDir 'modalSynthesisDiagramMSOutput'], '-depsc', '-r0')
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% LOOPBACK FM MS %%% 
    
    % loopback FM oscillators
    M = 175; % end index
    figure

    for i=1:Nf
        subplot(Nf, 1, i)
        plot(nVec(1:M), real(steelpan_z0_Mat(i,1:M)), 'linewidth', 2);
        xlim([nVec(1) nVec(M)]);
        set(gca,'linewidth', 3)
        set(gca,'XTick',[], 'YTick', [])
    end

    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    print([figDir 'modalSynthesisDiagramLoopbackFM'], '-depsc', '-r0')
    
    % Loopback FM output
    figure
    plot(nVec, real(steelpan_z0)/max(abs(real(steelpan_z0))), 'linewidth', 2);
    xlim([nVec(1) nVec(end)]);
    ylim([-1 1]);
    set(gca,'linewidth', 3)
    set(gca,'XTick',[], 'YTick', [])
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 1.6];
    print([figDir 'modalSynthesisDiagramLoopbackFMOutput'], '-depsc', '-r0')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% EXPONENTIALLY DECAYING AMPLITUDE ENVELOPES %%%
    for i=1:Nf
        subplot(Nf, 1, i)
        plot(nVec, env(i,:), 'linewidth', 2);
        xlim([nVec(1) nVec(end)]);
        ylim([0 1]);
        set(gca,'linewidth', 3)
        set(gca,'XTick',[], 'YTick', [])
    end
    
    fig = gcf
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 9];
    print([figDir 'modalSynthesisDiagramEnvelopes'], '-depsc', '-r0')
    
end
