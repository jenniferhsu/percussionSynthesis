This folder contains all the ways that I have come up with to produce interesting percussive sounds using my method.  

1. Run proofOfConcept.m first to generate modal synthesis percussion signals using the 3x3 mesh modal frequencies as well as theoretical steel plate modal frequencies.  We generate extended percussion signals using:
	a) rotational feedback FM
	b) stretched allpass filter formulated feedback FM
	c) time-varying allpass filters
	d) stretched allpass filter formulated feedback FM and time-varying allpass filter combination
	These sound examples are all written to a directory within audioExamples.  

2. Modify and run the percSynthTests.m file. This file uses convolutional synthesis to create more interesting percussive sounds.  You will need to modify the modal synthesis percussion signal location and the body resonator impulse response location within percSynthTests.m.  These files will be saved to audioExamples/convolutionalSynth/ directory

3. If you would like to run convolutional synthesis in a batch setting, use percSynthBatch.m and change the inputDir variable to your desired input directory.
