function yscaled = scaleForSavingAudio(y)
%SCALEFORSAVINGAUDIO scales audio input y for saving as a wav file without
% clipping. The maximum value will be +/- 0.95

yscaled = y/max(abs(y));
yscaled = 0.95 * yscaled;

end

