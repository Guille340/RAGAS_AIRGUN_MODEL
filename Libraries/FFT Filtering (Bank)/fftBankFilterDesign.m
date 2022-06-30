%  FftFilterBank = FFTBANKFILTERDESIGN(fs,bpo,freqLimits)
%
%  DESCRIPTION
%  Designs a bank of digital FFT filters to work with FFTBANKFILTER. The 
%  function returns a structure FFTFILTERBANK containing the sampling rate, 
%  cutoff frequencies and normalised cutoff frequencies of all filters in
%  the bank.
%
%  FFTFILTERBANK is used as input of FFTBANKFILTER to extract the metrics 
%  (RMS, exposure) from the bank filtered version of a  given broadband 
%  signal. For details about the content of the structure check the "OUTPUT 
%  VARIABLES" section.
%
%  FFTBANKFILTERDESIGN is a lightweight function used to perform error control 
%  over the design parameters of the FFT filter bank. This is done to minimise 
%  the number of operations carried out within FFTBANKFILTER, since the latter 
%  will be applied multiple times whereas the filter design will only happen 
%  once.
%
%  INPUT ARGUMENTS 
%  - fs: sampling rate on which the filter is based. Note that FS needs 
%    to be identical to the sampling rate of the signal to be filtered for 
%    cuttoff frequencies to be correct. Note that the signal to be filtered 
%    must have the same sampling rate if you want the filter to be applied 
%    at the correct frequencies.
%  - bpo: bandwidth factor. Positive integer that specifies the number of 
%    bands in an octave (i.e. BPO = 1, octave; BPO = 2, half-octave; BPO = 3, 
%    third-octave, etc). BPO determines the nominal bandwidth of the band-
%    pass filter. BPO is reciprocal of the bandwidth designator (see "Terms 
%    and Definitions" in ANSI S1.11 or BS EN 61260-1:2014).
%  - freqLimits: two-element vector indicating the bottom and top frequency 
%    limits of the filter bank [Hz].
%    
%  OUTPUT ARGUMENTS
%  - FftFilterBank: structure containing information for the FFT digital 
%    filter bank. FFTFILTERBANK comprises the following fields:
%    ¬ 'sampleRate': sampling rate [Hz] (see input variable FS). 
%    ¬ 'bandsPerOctave': bandwidth factor (see BPO input).
%    ¬ 'nominalFreq': nominal central frequencies of bands [Hz]
%    ¬ 'centralFreq': central frequencies of bands [Hz]
%    ¬ 'halfPowerFreq1': bottom cuttoff frequencies (-3dB) of bands [Hz]
%    ¬ 'halfPowerFreq2': top cutoff frequency (-3dB) [Hz]
%    ¬ 'centralFreqn': normalised central frequencies of bands [Hz]
%    ¬ 'halfPowerFreqn1': normalised bottom cutoff frequencies (-3dB) of 
%      bands [Hz]
%    ¬ 'halfPowerFreqn2': normalised top cutoff frequencies (-3dB) of
%      bands [Hz]
%
%  FUNCTION CALL
%  FftFilterBank = fftFilterBankDesign(fs,bpo,freqLimits) % ZEROPAD = TRUE
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
%  freqLimits = [20 10e3];
%
%  % 2) Design FFT Filter
%  FftFilterBank = fftBankFilterDesign(fs,bpo,freqLimits)
%
%  REFERENCES
%  - ANSI (2004), "ANSI S1.11: Specification for Octave, Half-Octave and
%    Third Octave Band Filter Sets", American National Standards Institute
%  - BSI (2014), "Electroacoustics - Octave-band and fractional-octave-band
%    filters. Part 1: Specifications". European Standard BS EN 61260-1:2014
%    published by the British Standards Institution (BSI).
%
%  See also FFTBANKFILTER

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  07 Aug 2021

function FftFilterBank = fftBankFilterDesign(fs,bpo,freqLimits)

% Check Number of Input Arguments
narginchk(3,3)

% Error Control (ARGUMENT: 'fs')
if numel(fs) > 1 || ~isreal(fs) || ~isscalar(fs) || fs < 0
    error('FS must be a real positive scalar number')
end

% Error Control (ARGUMENT: 'bpo')
if ~isnumeric(bpo) || ~isreal(bpo) || ~isscalar(bpo) || bpo < 1
    bpo = 1;
    warning(['BPO must be a real non-decimal scalar number equal to 1 '...
        'or higher. BPO = 1 will be used'])
elseif rem(bpo,1)
    bpo = round(bpo);
    warning('B must be a non-decimal number. BPO = %d will be used',bpo)
end

% Error Control (ARGUMENT: 'freqLimits')
if numel(freqLimits) == 2 && isreal(freqLimits)
    f1 = min(freqLimits); % bottom cutoff frequency [Hz]
    f2 = max(freqLimits); % top cutoff frequency [Hz]
    if isequal(f1,f2)
        f1 = 0;
        f2 = fs;
        warning(['FREQLIMITS must contain two different values. '...
            'FREQLIMITS set to [0 FS]'])
    end
    if f1 > f2
        fmin = f2;
        fmax = f1;
        f1 = fmin;
        f2 = fmax;
        warning(['The bottom cutoff frequency must be lower than the '...
            'top cutoff frequency (FREQLIMITS(1) < FREQLIMITS(2)). '...
            'The values will be swapped'])
    end
    if f1 < 0
        f1 = 0;
        warning('The frequency must be positive. FREQLIMITS(1) set to 0')
    end
    if f2 > fs/2
        f2 = fs/2;
        warning(['The top cutoff frequency must be lower than or equal to '...
            'half of the sample rate (FREQLIMITS(2) <= FS/2). '...
            'FREQLIMITS(2) set to FS/2'])
    end   
else
    error('FREQLIMITS must be a 2-element real vector')
end

% Calculate Cutoff Frequencies of Bandpass Filters
[fb,fbc_nom] = fractionalOctaveBands('BandsPerOctave',bpo,'FrequencyLimits',...
    [f1 f2],'LimitMode','bandedge','Base10',false); % central and edge frequencies of octave bands
fb1 = fb(:,1); % low-edge frequencies of fractional octave bands
fbc = fb(:,2); % central frequencies of fractional octave bands
fb2 = fb(:,3); % high-edge frequencies of fractional octave bands

% Filter Signal
FftFilterBank.sampleRate = fs;
FftFilterBank.bandsPerOctave = bpo;
FftFilterBank.nominalFreq = fbc_nom';
FftFilterBank.centralFreq = fbc';
FftFilterBank.halfPowerFreq1 = fb1';
FftFilterBank.halfPowerFreq2 = fb2';
FftFilterBank.centralFreqn = 2*fbc'/fs;
FftFilterBank.halfPowerFreqn1 = 2*fb1'/fs;
FftFilterBank.halfPowerFreqn2 = 2*fb2'/fs;
