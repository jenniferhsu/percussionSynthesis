This folder contains all the ways that I have come up with to produce interesting percussive sounds using my method.  

1. Run proofOfConcept.m first to generate modal synthesis percussion signals.  These percussion signals are modified with feedback FM and time-varying allpass filters and saved to the audio audioExamples/ directory.  

2. Modify and run the percSynthTests.m file. This file uses convolutional synthesis to create more interesting percussive sounds.  You will need to modify the modal synthesis percussion signal location and the body resonator impulse response location within percSynthTests.m.  These files will be saved to audioExamples/convolutionalSynth/ directory
