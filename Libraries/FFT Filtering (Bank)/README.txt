FFT Filtering (Bank)
=====================

MATLAB code for designing and applying a bank of bandpass filters in the
frequency domain using the Fast Fourier Transform (FFT).

The filtering function returns the metrics per band (SPLrms, SEL). 

Zero-padding is applied when the duration of the signal is comparable to or 
smaller than the period of the lowest frequency to be filtered, to avoid 
issues associated with low spectral resolution.

[Guillermo Jiménez Arranz, 12 May 2022]





