Digital Filtering (Single)
===========================

MATLAB code for designing and applying a single-band filter.

The filtering function can either return the filtered signal or the 
corresponding metrics (SPLrms, SEL, SPLpk, SPLpp). 

Zero-padding is applied when the duration of the signal is comparable to or 
smaller than the period of the lowest frequency to be filtered, to avoid 
issues associated with delay of filters. 

Zero-phase filtering can be applied to keep the phase information of the 
signal intact. This can be useful for estimating the peak amplitude of 
filtered bands.

[Guillermo Jiménez Arranz, 12 May 2022]





