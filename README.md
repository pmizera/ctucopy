# CtuCopy 
# Speech Enhancement, Feature Extraction tool

CtuCopy is a command line tool implementing speech enhancement and feature extraction algorithms. It is similar to HCopy tool from Cambridge's HTK and is partially compatible with it. CtuCopy acts as a filter with speech waveform file(s) at the input and either a speech waveform file(s) or feature file(s) at the output.

CtuCopy implements several speech enhancing methods based on spectral subtraction. Extended Spectral Subtraction - exten [Sovka et al., EUSIPCO'96] combines Wiener filtering and spectral subtraction with no need of a voice activity detector (VAD). Other methods are based on spectral subtraction with VAD which can be either external (data read from file) or internal (Burg's cepstral detector). The noise suppression methods can be applied either directly to speech spectra or to a spectra filtered by a bank of filters.

Bank of filters offers a set of frequency scales (melodic, Bark, expolog, linear) and various filter shapes (triangular, rectangular, trapezoidal same as time-domain IIR based implementation) which can be combined to form standard and user-defined filter banks.

In feature extraction mode a number of common features can be extracted from original or enhanced speech, e.g. LPC, PLP, MFCC or magnitude spectrum. Current version also implements TRAP-DCT features. Features can be saved either in HTK, pfile, or KALDI format.
