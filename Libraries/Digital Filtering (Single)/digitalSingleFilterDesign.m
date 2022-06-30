%  DigitalFilter = DIGITALSINGLEFILTERDESIGN(fs,freqCutoff,varargin)
%
%  DESCRIPTION
%  Designs a digital zero-phase IIR filter to work with DIGITALSINGLEFILTER. 
%  The function returns a structure DIGITALFILTER containing the filtering 
%  object and information about the decimator and digital filter. 
%  DIGITALFILTER is used as input of DIGITALSINGLEFILTER to extract the 
%  filtered signals or metrics (RMS, exposure, peak, peak-to-peak) from a 
%  given broadband signal. For details about the content of the structure 
%  check the "OUTPUT VARIABLES" section.
% 
%  INPUT ARGUMENTS (Fixed)
%  - fs: sampling rate on which the filter is based. Note that FS needs 
%    to be identical to the sampling rate of the signal to be filtered for 
%    cuttoff frequencies to be correct. Note that the signal to be filtered 
%    must have the same sampling rate if we want the filter to be applied at 
%    the correct frequencies.
%  - freqCutoff: two-element vector containing the bottom and top cutoff 
%    frequencies [Hz]. The top frequency must be lower or equal than 
%    the Nyquist frequency (FREQLIMITSS(2) <= FS/2. The bottom frequency
%    must be higher or equal than zero (FREQCUTOFF(1) >= 0). FREQCUTOFF
%    are additionally used to identify the filter type: for a low-pass filter
%    (LPF), set FREQCUTOFF(1) = 0; for a high-pass filter (HPF), set
%    FREQCUTOFF(2) = FS/2; for a band-pass filter (BPF), set FREQCUTOFF(1) > 0 
%    and FREQCUTOFF(2) < FS/2.
%
%  INPUT ARGUMENTS (Variables, Property/Value Pairs)
%  In the function call, type the property string followed by its value (comma-
%  separated). Property/value pairs are variable input arguments specified after
%  the fixed arguments. Any number of property/value pairs can be included in 
%  the function call. The following propeties are available.
%  - 'FilterOrder': order of the decimator and bandpass filters. Both
%     types of filters use an IIR Unconstrained Butterworth response. The 
%     order must be even and >= 2. Use FILTORDER >= 4 for reliable bandpass 
%     filtering. If this property is omitted, the default FILTORDER = 8 is
%     used.
%  - 'FilterType': type of IIR filter. Three options: 'lowpass','highpass',
%     'bandpass' and 'argmag'. A fifth option ('bandstop') may be included in 
%     a future version. If FILTERTYPE option is omitted, 'bandpass' will be 
%     assumed.
%  - 'Amplitudes': numeric vector representing the magnitude response of the
%     arbitrary magnitude filter 'arbmag'.
%   
%  OUTPUT VARIABLES
%  - DigitalFilter: structure containing information and filtering objects 
%    for the decimator and target filter. DIGITALFILTER comprises the 
%    following fields and subfields:
%
%    ~~~ Level 1 Fields (DigitalFilter.<Level1Field>) ~~~
%    ¬ 'sampleRate': sampling rate [Hz] (see input variable FS). 
%    ¬ 'freqCutoff': two-element vector indicating the bottom and top frequency 
%      limits of the filter bank [Hz] (see input variable FREQCUTOFF).
%    ¬ decimator: substructure containing the filtering object and specific
%      information for the unique decimator used in the filter bank design.
%    ¬ target: substructure containing the filtering object designed
%      with MATLAB's DESIGNFILT.
%
%    ~~~ Level 2 Fields (decimator.<Level2Field>) ~~~
%    ¬ 'filterType': 'lowpass' for decimation filter.
%    ¬ 'filterResponse': 'IIR' or infinite impulse response (no FIR filters)
%    ¬ 'filterOrder': order of the decimator. The filter order is given 
%      specifically as an input of DIGITALSINGLEFILTERDESIGN or calculated 
%      from the target filter class. The same order is used for all filters 
%      in the bank, including the decimator.
%    ¬ 'filterObject': digitalFilter class for an unconstrained, lowpass, IIR
%      Butterworth filter. This is the anti-aliasing filter to be applied
%      before downsampling. The filter is designed for its bandedge frequency
%      to be 1.5 times the central frequency of the octave band to be filtered.
%      Decimation by 2 is applied when the central frequency of the octave band
%      to be filtered is at least 2^2.5 times lower than the current sampling 
%      rate (i.e. sampling rate before downsampling by 2). This means that the 
%      normalised frequencies of the filter (and the filter itself) are the 
%      same for all decimation stages, therefore only one decimation filter is 
%      needed. A signal can be filtered with this object using the functions 
%      FILTER or FILTFILT (see DIGITALSINGLEFILTER).
%    ¬ 'halfPowerFreqn': normalised half-power (-3 dB) frequency for the
%      decimation filter. The absolute frequency is calculated by
%      multiplying HALFPOPWERFREQN * FR/2, where FR is the sampling rate
%      for a given decimation stage (FR = 2*FS/D, before downsampling by 2).
%    ¬ 'decimationFactor': vector of factors by which the original sampling 
%      frequency FS and the number of samples in the original signal are 
%      reduced at each decimation stage (after downsampling by 2). As many 
%      elements as octave bands containing at least one fractional octave band.
%    ¬ 'sampleRate': sampling rate at each decimation stage after downsampling 
%      by 2 (FR = FS/D).
%    ¬ 'halfPowerFreq': absolute half-power (-3 dB) frequency for the
%      decimation filter at each decimation stage. As many elements as octave
%      bands containing at least one fractional octave band. HALFPOWERFREQ
%      is NaN for the parent octave bands that do not have a decimation filter.
%    
%    ~~~ Level 2 Fields (target.<Level2Field>) ~~~
%    ¬ filterType: 'bandpass','lowpass',or 'highpass'.
%    ¬ filterResponse: 'IIR' or infinite impulse response (no FIR filters)
%    ¬ filterOrder: order of the unique bandpass filters. The filter order is 
%      given specifically as an input of DIGITALSINGLEFILTERDESIGN.
%    ¬ filterObject: unconstrained IIR Butterworth filter of digitalFilter 
%      class. This filter corresponds to the filter to be applied after
%      decimation. A signal can be filtered with this filter object using the 
%      functions FILTER or FILTFILT (see DIGITALSINGLEFILTER).
%    ¬ halfPowerFreqn1: normalised low-edge frequency. The absolute frequency
%      is calculated by multiplying CUTOFFFREQN1 * FR/2, where FR is the 
%      sampling rate for the last decimation stage (FR = FS/D, after 
%      downsampling by 2).
%    ¬ halfPowerFreqn2: normalised high-edge frequency. he absolute frequency 
%      is calculated by multiplying CUTOFFFREQN1 * FR/2, where FR is the 
%      sampling rate for the last decimation stage (FR = FS/D, after 
%      downsampling by 2).
%    ¬ sampleRate: sample rate at the last decimation stage after downsampling 
%      by D = 2^N, with N the number of decimation stages (FR = FS/D).
%    ¬ halfPowerFreq1: absolute low-edge frequency.
%    ¬ halfPowerFreq2: absolute high-edge frequency.
%
%  CONSIDERATIONS & LIMITATIONS
%  - The FILTERTYPE property is currently redundant, since the type of
%    filter can be directly inferred from FREQUENCYLIMITS. However, the
%    property is left as anticipation for a future revision that may
%    include band-stop filtering, in which case the FILTERTYPE property
%    will be necessary to differentiate between 'bandpass' and 'bandstop'.
%
%  FUNCTION CALL
%  1. DigitalFilter = DIGITALSINGLEFILTERDESIGN(fs,freqCutoff)
%  2. DigitalFilter = DIGITALSINGLEFILTERDESIGN(...,PROPERTYNAME,PROPERTYVALUE)
%
%  If no properties are specified, their default values will be used. The 
%  default values are: FILTERORDER = 8; FILTERTYPE = 'bandpass'; AMPLITUDES
%  = [1 1].
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
%  filtOrder = 10;
%
%  % 2) Filter Bank Design
%  DigitalFilter = DIGITALSINGLEFILTERDESIGN(fs,freqCutoff,'FilterOrder',...
%      filtOrder)
%
%  See also DIGITALSINGLEFILTER

%  VERSION 2.1
%  Date: 04 Mar 2022
%  Author: Guillermo Jimenez Arranz
%  - Added 'arbmag' filter to FILTERTYPE property.
%  - Added 'amplitude' property for 'arbmag' filter.
%
%  VERSION 2.0
%  Date: 06 Aug 2021
%  Author: Guillermo Jimenez Arranz
%  - Renamed functions to inclulde the word 'single'. This is done to avoid
%    conflict with the ANSI bank filtering toolbox.
%
%  VERSION 1.2
%  Date: 21 May 2021
%  Author: Guillermo Jimenez Arranz
%  - Updated help
%  - Removed 'bandstop' option in FILTERTYPE property as it cannot be 
%    easily implemented. In a future revision, DIGITALSINGLEFILTERDESIGN could include
%    a 'bandstop' option, where the band stop filter is implemented as a
%    highpass and lowpass filter.
%
%  VERSION 1.1
%  Date: 16 May 2021
%  Author: Guillermo Jimenez Arranz
%  - Corrected the normalised half-power frequency of the decimator, from
%    2*fd1/fr to fd1/fr (line 264).
%  - Added variable input argument 'FILTERTYPE'.
%  - Added 'BandStop' filter response type to the available 'LowPass',
%    'HighPass', and 'Bandpass'.
%  - Updated the Error Control section.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  26 Mar 2021

function DigitalFilter = digitalSingleFilterDesign(fs,freqCutoff,varargin)

% Check Number of Input Arguments
nFixArg = 2;
nProArg = 6;
narginchk(nFixArg,nFixArg + nProArg)
if rem(nargin - nFixArg,2)
    error('Property and value input arguments must come in pairs')
end

% Extract and Verify Input Properties
properties_valid = {'filterorder','filtertype','amplitudes'};
properties = lower(varargin(1:2:end));
if any(~ismember(properties,properties_valid))
    error('One or more input properties are not recognised')
end

% Default Input Values
filtOrder = 8; % filter order
filtType = 'bandpass'; % filter type
amp = [1 1]; % amplitudes of arbitrary filter (flat response)

% Extract and Verify Input Values
values = varargin(2:2:end);
nPairs = (nargin - nFixArg)/2; % number of (PROPERTY,VALUE) pairs
for m = 1:nPairs
    property = properties{m};
    switch property % add new properties with CASE PROPERTY
        case 'filterorder'
            filtOrder = values{m};
        case 'filtertype'
            filtType = values{m};
        case 'amplitudes'
            amp = values{m};
    end
end

% Set Cutoff Frequencies to 0 and FS/2 for FILTTYPE = 'arbmag'
if strcmp(filtType,'arbmag')
    freqCutoff(1) = 0;
    freqCutoff(2) = fs/2;
end

% Error Control: FS
if numel(fs) > 1 || ~isreal(fs) || fs < 0
    error('FS must be a real positive number')
end

% Error Control: FILTTYPE
if ~ismember(filtType,{'lowpass','highpass','bandpass','arbmag'})
    error('Invalid argument for ''FilterType'' input property')
end

% Error Control: FREQCUTOFF
if numel(freqCutoff) == 2 && isreal(freqCutoff)
    f1 = min(freqCutoff); % bottom cutoff frequency [Hz]
    f2 = max(freqCutoff); % top cutoff frequency [Hz]
    if isequal(f1,f2)
        f1 = 0;
        f2 = fs/2;
        warning(['FREQCUTOFF must contain two different values. '...
            'FREQCUTOFF set to [0 FS/2]'])
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
        warning('The frequency  must be positive. FREQCUTOFF(1) set to 0')
    end
    if f2 > fs/2
        f2 = fs/2;
        warning(['The top cutoff frequency must be lower than or equal '...
            'to half of the sample rate (FREQCUTOFF(2) <= FS/2). '...
            'FREQCUTOFF(2) set to FS/2'])
    end   
    if f1 == 0 && ~strcmp(filtType,'lowpass') % bottom cutoff frequency < 0 
        if f2 < fs/2
            filtType = 'lowpass';
            warning(['MIN(FREQCUTOFF) = 0 and MAX(FREQCUTOFF) < FS/2. '...
                'A low-pass filter will be applied.'])
        else % f2 == fs/2
            warning(['MIN(FREQCUTOFF) = 0 and MAX(FREQCUTOFF) = FS/2. '...
                'No filtering will be applied'])
        end    
    end
    if f2 == fs/2 && ~strcmp(filtType,'highpass') % bottom cutoff frequency < 0 
        if f1 > 0
            filtType = 'highpass';
            warning(['MIN(FREQCUTOFF) > 0 and MAX(FREQCUTOFF) == FS/2. '...
            'A high-pass filter will be applied.'])
        else % f1 == 0
            warning(['MIN(FREQCUTOFF) = 0 and MAX(FREQCUTOFF) = FS/2. '...
                'No filtering will be applied'])
        end    
    end
    if f1 > 0 && f2 < fs/2 && ~ismember(filtType,{'bandpass'})
        filtType = 'bandpass';
        warning(['MIN(FREQCUTOFF) > 0 and MAX(FREQCUTOFF) < FS/2. '
        'A band-pass filter will be applied.'])
    end
else
    error('FREQCUTOFF must be a 2-element real vector')
end

% Error Control: AMP
if strcmp(filtType,'arbmag') && (~isnumeric(amp) || ~isvector(amp))
    amp = [1 1];
    warning(['Non-supported value for input PROPERTY = ''Amplitudes''.'...
        'No filtering will be applied']) 
end

% Error Control: FILTORDER
if ~isnumeric(filtOrder) || ~isscalar(filtOrder) || ~isreal(filtOrder)...
        || rem(filtOrder,2) || filtOrder < 2 
    filtOrder = round(filtOrder/2)*2; % even filter order
    if filtOrder < 2
        filtOrder = 8;
        warning(['Input property ''FiltOrder'' must be an even number '...
            'greater than or equal to 2. FILTORDER = %d will be used'],...
            filtOrder)
    end
end
if strcmp(filtType,'arbmag') && filtOrder < 2*length(amp)
    warning(['FILTORDER >= 2*LENGTH(AMPLITUDES) is recommended for an '...
        'arbitrary magnitude filter (FILTERTYPE = ''arbmag''). Check with '...
        'FVTOOL that the filter response is close to the specified '...
        'magnitude response'])
end

% Decimation Factor and Sample Rate after Decimation
fd1 = f2 * 2; % half-power frequency of current decimator
frMin = min([fs, 2*fd1]); % minimum sampling rate to avoid aliasing
D = 2^floor(log2(fs/frMin)); % decimation factor for current octave (after downsampling by 2)
fr = fs/D; % decimated sampling rate (after downsampling by 2)
if D == 1, fd1 = NaN; end % no decimation filter for unity decimation factor

% Initialise Parameters of Decimation Filter
decimFiltType = 'lowpass'; % type of filter
decimFiltResponse = 'IIR'; % impulse response of filter
decimFiltOrder = 0; % filter order
decimFilter = []; % filter object
decimFn1 = []; % normalised half-power frequency
decimDecimFactor = []; % decimation factor (after downsampling by 2)
decimFd = fs; % sampling rate (before downsampling by 2)
decimF1 = []; % half-power frequencies

% Parameters of Decimator
if ~isnan(fd1)
    decimFiltOrder = filtOrder;
    decimFn1 = fd1/fr; % normalised half-power frequency (before downsampling by 2)
    decimFilter = designfilt('lowpassiir','FilterOrder',decimFiltOrder,...
                'HalfPowerFrequency',decimFn1); % design decimator (anti-alias filter)
    decimDecimFactor = D; % decimation factor (after downsampling by 2)
    decimFd = min([fs, fr]); % sampling rate (after downsampling by 2)
    decimF1 = fd1; % absolute half-power frequency of decimator
end

% Initialise Parameters of Bandpass Filter
targetFiltType = ''; % type of filter
targetFiltResponse = ''; % impulse response of filter
targetFiltOrder = 0; % filter order
targetFilter = []; % filter object
targetFn1 = []; % normalised low-edge half-power frequency
targetFn2 = []; % normalised high-edge half-power frequency
targetFd = fs; % sampling rate (after downsampling by 2)
targetF1 = []; % low-edge half-power frequencies
targetF2 = []; % high-edge half-power frequencies

% Design Target Filter
if f1 ~= 0 || f2 ~= fs/2 % filtType = {'lowpass','highpass','bandpass'}
    targetFiltType = filtType; % type of filter
    targetFiltResponse = 'IIR'; % impulse response of filter
    targetFn1 = 2 * f1/fr; % normalised low-edge frequency of current bandpass filter
    targetFn2 = 2 * f2/fr; % normalised high-edge frequency of current bandpass filter
    targetFiltOrder = filtOrder;
    targetFd = fr; % sampling rates (after downsampling by 2)
    targetF1 = f1; % low-edge frequency for bandpass filters
    targetF2 = f2; % high-edge frequency for bandpass filters
    switch filtType
        case 'lowpass'
            targetFilter = designfilt('lowpassiir','FilterOrder',...
                targetFiltOrder,'HalfPowerFrequency',targetFn2);
        case 'highpass'
            targetFilter = designfilt('highpassiir','FilterOrder',...
                targetFiltOrder,'HalfPowerFrequency',targetFn1); 
        case 'bandpass'
            targetFilter = designfilt('bandpassiir','FilterOrder',...
                targetFiltOrder,'HalfPowerFrequency1',targetFn1,...
                'HalfPowerFrequency2',targetFn2); 
    end
    
elseif strcmp(filtType,'arbmag') % filtType = 'arbmag'
    targetFiltType = 'arbmag'; % type of filter
    targetFiltResponse = 'FIR'; % impulse response of filter
    targetFiltOrder = filtOrder; % filter order
    targetFn1 = 2 * f1/fs; % normalised low-edge half-power frequency
    targetFn2 = 2 * f2/fs; % normalised high-edge half-power frequency
    targetFd = fs; % sampling rate (after downsampling by 2)
    targetF1 = f1; % low-edge half-power frequencies
    targetF2 = f2; % high-edge half-power frequencies
    targetFilter = designfilt('arbmagfir','FilterOrder',targetFiltOrder,...
            'Frequencies',linspace(0,1,length(amp)),'Amplitudes',amp);  
end

% Store General Parameters in ¦filterStructure¦
DigitalFilter.sampleRate = fs; % original sampling frequency [Hz]
DigitalFilter.freqCutoff = freqCutoff; % frequency range [fmin fmax]

% Store Decimator Filterbank in ¦filterStructure¦
DigitalFilter.decimator.filterType = decimFiltType;
DigitalFilter.decimator.filterResponse = decimFiltResponse;
DigitalFilter.decimator.filterOrder = decimFiltOrder;
DigitalFilter.decimator.filterObject = decimFilter;
DigitalFilter.decimator.halfPowerFreqn = decimFn1;
DigitalFilter.decimator.decimationFactor = decimDecimFactor;
DigitalFilter.decimator.sampleRate = decimFd;
DigitalFilter.decimator.halfPowerFreq = decimF1;

% Store Bandpass Filterbank in ¦filterStructure¦
DigitalFilter.target.filterType = targetFiltType;
DigitalFilter.target.filterResponse = targetFiltResponse;
DigitalFilter.target.filterOrder = targetFiltOrder;
DigitalFilter.target.filterObject = targetFilter;
DigitalFilter.target.halfPowerFreqn1 = targetFn1;
DigitalFilter.target.halfPowerFreqn2 = targetFn2;
DigitalFilter.target.sampleRate = targetFd;
DigitalFilter.target.halfPowerFreq1 = targetF1;
DigitalFilter.target.halfPowerFreq2 = targetF2;
