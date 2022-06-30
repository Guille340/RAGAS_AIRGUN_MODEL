%  FftFilter = FFTSINGLEFILTERDESIGN(fs,freqCutoff)
%
%  DESCRIPTION
%  Designs a digital FFT filter to work with FFTSINGLEFILTER. The function 
%  returns a structure FFTFILTER containing the sampling rate, cutoff
%  frequencies and normalised cutoff frequencies.
%
%  FFTFILTER is used as input of FFTSINGLEFILTER to extract the metrics 
%  (RMS, exposure) from the filtered version of a given broadband signal. 
%  For details about the content of the structure check the "OUTPUT ARGUMENTS" 
%  section.
%
%  FFTSINGLEFILTERDESIGN is a lightweight function used to perform error 
%  control over the design parameters of the FFT filter. This is done to 
%  minimise the number of operations carried out within FFTSINGLEFILTER, 
%  since the latter will be applied multiple times whereas the filter design
%  will only happen once.
%
%  INPUT ARGUMENTS 
%  - fs: sampling rate on which the filter is based [Hz]. Note that FS needs 
%    to be identical to the sampling rate of the signal to be filtered for 
%    cuttoff frequencies to be correct.
%  - freqCutoff: two-element vector containing the bottom and top cutoff 
%    frequencies of the filter [Hz]. The top frequency must be lower than 
%    or equal to the Nyquist frequency (FREQLIMITSS(2) <= FS/2. The bottom 
%    frequency must be higher than or equal to zero (FREQCUTOFF(1) >= 0).
%    
%  OUTPUT ARGUMENTS
%  - FftFilter: structure containing information for the FFT digital filter.
%    FFTFILTER comprises the following fields:
%    ¬ 'sampleRate': sampling rate used to design the filter [Hz] (see input 
%      variable FS). Note that the signal to be filtered must have the 
%      same sampling rate if we want the filter to be applied at the correct 
%      frequencies.
%    ¬ 'halfPowerFreq1': bottom cuttoff frequency (-3dB) [Hz]
%    ¬ 'halfPowerFreq2': top cutoff frequency (-3dB) [Hz]
%    ¬ 'halfPowerFreqn1': normalised bottom cutoff frequency (-3dB)
%    ¬ 'halfPowerFreqn2': normalised top cutoff frequency (-3dB)
%
%  FUNCTION CALL
%  FftFilter = FFTSINGLEFILTERDESIGN(fs,freqCutoff) % ZEROPAD = TRUE
%
%  FUNCTION DEPENDENCIES
%  - None
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  EXAMPLE
%  % 1) Configuration Data
%  fs = 44100;
%  freqCutoff = [20 10e3];
%
%  % 2) Design FFT Filter
%  FftFilter = fftSingleFilterDesign(fs,freqCutoff)
%
%  See also FFTSINGLEFILTER

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  07 Aug 2021

function FftFilter = fftSingleFilterDesign(fs,freqCutoff)

% Check Number of Input Arguments
narginchk(2,2)

% Error Control (ARGUMENT: 'fs')
if numel(fs) > 1 || ~isreal(fs) || ~isscalar(fs) || fs < 0
    error('FS must be a real positive scalar number')
end

% Error Control (ARGUMENT: 'freqCutoff')
if numel(freqCutoff) == 2 && isreal(freqCutoff)
    f1 = min(freqCutoff); % bottom cutoff frequency [Hz]
    f2 = max(freqCutoff); % top cutoff frequency [Hz]
    if isequal(f1,f2)
        f1 = 0;
        f2 = fs;
        warning(['FREQCUTOFF must contain two different values. '...
            'FREQCUTOFF set to [0 FS]'])
    end
    if f1 > f2
        fmin = f2;
        fmax = f1;
        f1 = fmin;
        f2 = fmax;
        warning(['The bottom cutoff frequency must be lower than the '...
            'top cutoff frequency (FREQCUTOFF(1) < FREQCUTOFF(2)). '...
            'The values will be swapped'])
    end
    if f1 < 0
        f1 = 0;
        warning('The frequency must be positive. FREQCUTOFF(1) set to 0')
    end
    if f2 > fs/2
        f2 = fs/2;
        warning(['The top cutoff frequency must be lower than or equal to '...
            'half of the sample rate (FREQCUTOFF(2) <= FS/2). '...
            'FREQCUTOFF(2) set to FS/2'])
    end   
else
    error('FREQCUTOFF must be a 2-element real vector')
end

% Filter Signal
FftFilter.sampleRate = fs;
FftFilter.halfPowerFreq1 = f1;
FftFilter.halfPowerFreq2 = f2;
FftFilter.halfPowerFreqn1 = 2*f1/fs;
FftFilter.halfPowerFreqn2 = 2*f2/fs;
