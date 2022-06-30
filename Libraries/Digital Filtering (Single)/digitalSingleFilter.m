%  xb = DIGITALSINGLEFILTER(DigitalFilter,x,varargin)
%
%  DESCRIPTION
%  Filters the input signal X using a filter defined in DIGITALSINGLEFILTER 
%  structure. DIGITALSINGLEFILTER applies an appropriate amount of zero-
%  padding to the signal. The output XB is a zero-phase filtered version of 
%  X. By using zero-padding, DIGITALSINGLEFILTER eliminates any issue that 
%  may arise when filtering signals with a length comparable to the filter’s 
%  group delay. 
%
%  INPUT ARGUMENTS
%  - DigitalFilter: filtering structure made with DIGITALSINGLEFILTERDESIGN. 
%    The structure contains information about the decimation and target 
%    filter, and their corresponding filter objects generated with DESIGNFILT. 
%  - x: vector of data to be filtered. Its length must be three times
%    larger, or more, than the filter order (if unknown, use FILTORD to 
%    extract the filter order).
%  - PROPERTIES (varargin): the user can input up to four function properties. 
%    The properties must be given in (Name,Value) pairs, separated by comma 
%    with the property name first. The available property names and the 
%    description of their values are given below:
%    ¬ 'MetricsOutput': logical or numeric [0,1] value indicating whether the 
%       output variable XB shall return the waveform or metrics (RMS, Exposure, 
%       Peak, Peak to Peak) or the filtered signal. If this property is 
%       omitted, the default METRICSOUTPUT = FALSE is used.
%       # FALSE (or 0): XB is a vector representing the waveform of the signal 
%         filtered with DIGITALSINGLEFILTER. The length of XB is that of the 
%         zero-padded signal if ZEROPADDED = TRUE and DATAWRAP = FALSE. 
%         Otherwise, the length of XB will be the same as X. Slower than 
%         METRICSOUTPUT = TRUE, but gives the filtered waveforms for further 
%         processing. DEFAULT option.
%       # TRUE (or 1): XB is a four-element vector containing the amplitude/
%         energy metrics of the filtered signal obtained with the filter 
%         designed with DESIGNFILT. The vector contains the following metrics, 
%         in ascending order: root-mean square (RMS), exposure, peak and 
%         peak-to-peak. Note that the units are those of the original signal.
%    ¬ 'ZeroPadding': logical or [0,1] numeric value indicating whether zero-
%       padding will be applied to the input signal X. Zero-padding is
%       necessary in those cases where the signal is short enough that lowpass
%       or bandpass filtering will displace (delay) the filtered signal beyond 
%       its original window length. This will result in a noticeable energy 
%       error. Depending on the signal duration and top cutoff frequency of 
%       the filter, the error can be as large as several tens of decibels. 
%       One-sided zero padding is applied when FILTERMODE = 'filter', and 
%       two-sided zero padding when using FILTERMODE = 'filtfilt'. Any 
%       filtering approach will suffer from this "time leakage" effect. If 
%       this property is omitted, ZEROPADDING = TRUE (default).
%     ¬ 'FilterMode': character string that specifies the function used for
%       filtering. If this property is omitted, the default FILTERMODE = 
%       'filter' is used. Two options:
%       # 'filter': fastest option (DEFAULT).
%       # 'filtfilt': zero-phase filtering option for IIR filters. Slower
%         than 'filter'. Use if the phase of the filtered signal is important 
%         and want its envelope to follow the envelope of the original signal 
%         (for example, for peak amplitude calculations). Note that custom 
%         function MYFILTFILT is used instead of MATLAB's FILTFILT due to large 
%         amplitude errors; these errors are caused by FILTFILT trying to 
%         minimise transient effects at the edges of the signal to be filtered.
%     ¬ 'DataWrap': logical or [0,1] numeric value indicating that wrapping
%        will be applied to the bandpass filtered signals in XB when choosing
%        the option 'METRICSOUTPUT' = FALSE. The function uses MATLAB's 
%        DATAWRAP to reduce the zero-padded signal ('ZEROPADDING' = TRUE) to 
%        the length of the original input signal X. The wrapping will result 
%        in a certain signal energy error due to amplitude cancellations, in 
%        particular at the lowest frequencies. Those minor errors are 
%        compensated by applying a correction factor to the wrapped signals.
%        If this property is omitted, DATAWRAP = TRUE (default). This property 
%        will be ignored when 'METRICSOUTPUT' = TRUE.
% 
%  OUTPUT ARGUMENTS
%  - y: zero-phase filtered output. Alternatively, it will contain the 
%    amplitudes for the selected metrics (srms, sexp, sz2p, sp2p).
%
%  FUNCTION CALL
%  xb = DIGITALSINGLEFILTER(DigitalFilter,x)
%  xb = DIGITALSINGLEFILTER(...,PROPERTYNAME,PROPERTYVALUE)
%
%  If no properties are specified, their default values will be used. The list 
%  of default values is:
%
%   METRICSOUTPUT = 0;
%   ZEROPADDING = TRUE;
%   FILTERMODE = 'filter';
%   DATAWRAP = TRUE;
%
%  FUNCTION DEPENDENCIES
%  - digitalSingleFilterDesign
%  - myfiltfilt
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  EXAMPLE
%  % 1) Configuration Data
%  fs = 44100;
%  freqLims = [20 10e3];
%  filtOrder = 10;
%  T = 1; % signal duration [s]
%  x = 2*rand(1,T*fs) - 1; % signal
%  filtMode = 'filtfilt';
%
%  % 2) Filter Bank Design
%  DigitalFilter = digitalSingleFilterDesign(fs,freqLims,'FilterOrder',...
%      filtOrder)
%
%  % 3) Filter Signal with IIR and FFT Methods
%  xb = digitalSingleFilter(DigitalFilter,x,'MetricsOutput',true,...
%      'FilterMode',filtMode,'DataWrap',true);
%
%  % 4) Band Levels
%  Xb_rms = 20*log10(xb(1));
%  Xb_exp = 10*log10(xb(2));
%  Xb_z2p = 20*log10(xb(3));
%  Xb_p2p = 20*log10(xb(4));
%
%  See also DIGITALSINGLEFILTERDESIGN, MYFILTFILT 

%  VERSION 2.1
%  Date: 04 Nov 2021
%  Revised by: Guillermo Jimenez Arranz
%  - Use of ISDIGITALSINGLEFILTER for checking the validity of the filter.
%
%  VERSION 2.0
%  Date: 06 Aug 2021
%  Revised by: Guillermo Jimenez Arranz
%  - Renamed functions to inclulde the word 'single'. This is done to avoid
%    conflict with the ANSI bank filtering toolbox.
%  - Removed the FFT filtering functionality. This will be implemented
%    separately for both bank and single filtering.
%
%  VERSION 1.2
%  Date: 21 May 2021
%  Revised by: Guillermo Jimenez Arranz
%  - Updated help
%  - Removed 'bandstop' option in FILTERTYPE property as it cannot be 
%    easily implemented. In a future revision, DIGITALSINGLEFILTER could 
%    include a 'bandstop' option, where the band stop filter is implemented 
%    as a highpass and lowpass filter.
%
%  VERSION 1.1
%  Date: 16 May 2021
%  Revised by: Guillermo Jimenez Arranz
%  - Added FFT filtering.
%  - Added 'DataWrap' input option to reduce the zero-padded signal to the
%    length of the original input signal X.
%  - Included SEL, peak, and peak-to-peak metrics in XB when 'MetricsOutput'
%    1 or 2 are selected.
%  - Replaced FILTFILT function with MYFILTFILT due to problems associated
%    with the former that result in large transients at the start and end
%    of the signal.
%  - Added output for empty (flat) filter.
%  - Corrected interpolation in METRICSOUTPUT = 0 to account for the
%    possibility of having to pad with zeroes (line 304-323).
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  27 Mar 2021

function xb = digitalSingleFilter(DigitalFilter,x,varargin)

% Check Number of Input Arguments
narginchk(2,10)
nVarargin = nargin - 2;
if rem(nVarargin,2)
    error('Property and value input arguments must come in pairs')
end

% Initialise Default Parameters
metricsOutput = 0; 
zeroPad = true;
dataWrap = true;
filtMode = 'filter';

% Retrieve Input Variables
for m = 1:2:nVarargin
    inputProperty = lower(varargin{m}); % case insensitive
    inputProperties = lower({'MetricsOutput','ZeroPadding','FilterMode',...
        'DataWrap'});
    if ~ismember(inputProperty,inputProperties)
        error('Invalid input property')
    else
        switch inputProperty
            case 'metricsoutput'
                metricsOutput = varargin{m+1};
            case 'zeropadding'
                zeroPad = varargin{m+1};
            case 'filtermode'
                filtMode = varargin{m+1};
            case 'datawrap'
                dataWrap = varargin{m+1};
        end
    end
end

% Error Control (check filter structure)
if ~isDigitalSingleFilter(DigitalFilter)
    error('The input argument DIGITALSINGLEFILTER is not a valid filter structure')
end

% Error Control (check input signal)
if ~isreal(x) || ~isvector(x)
    error('The input signal must be a real numeric vector')
end

% Error Control (PROPERTY: 'MetricsOutput')
if all(metricsOutput ~= [0 1 2])
    metricsOutput = 0;
    warning(['Non-valid value for METRICSOUTPUT property. '...
        'METRICSOUTPUT = 0 will be used'])
end

% Error Control (PROPERTY: 'ZeroPadding')
if ~islogical(zeroPad) || all(zeroPad ~= [0 1])
    zeroPad = true;
    warning(['ZEROPADDING property must be [0,1] or logical. '...
        'ZEROPADDING = TRUE will be used'])
end

% Error Control (PROPERTY: 'FilterMode')
switch filtMode
    case 'filter'
        filtFun = @filter;
    case 'filtfilt'
        filtFun = @myfiltfilt;
    otherwise
        filtFun = @filter;
        warning(['Non-valid value for FILTERMODE property. FILTERMODE = '...
            '''filter'' will be used'])
end

% Error Control (PROPERTY: 'DataWrap')
if ~islogical(dataWrap) || all(dataWrap ~= [0 1])
    dataWrap = true;
    warning(['DATAWRAP property must be [0,1] or logical. '...
        'DATAWRAP = TRUE will be used'])
end

% Load Parameters
fs = DigitalFilter.sampleRate; % sampling rate of filter and input signal
decimField = DigitalFilter.decimator; % substructure of decimation filter
targetField = DigitalFilter.target; % substructure of target filter
f1 = targetField.halfPowerFreq1; % low cutoff frequency of target filter
f2 = targetField.halfPowerFreq2; % high cutoff frequency of target filter
targetFilter = targetField.filterObject; % target filter object
targetFilterType = targetField.filterType; % type of target filter ('lowpass','highpass','bandpass')
decimFilter = decimField.filterObject; % decimation filter object
decimFactor = decimField.decimationFactor; % decimation factor
        
% Filter Input Signal
x = x(:); % convert x into column vector
L = length(x); % length of input signal [samples]
Lz = L; % intialise signal length to no zero-padding [samples]
if ~isempty(targetFilterType) % filter if filter exists
    if metricsOutput < 2 % Filtering with Filter
        % Zero-Padding
        if zeroPad        
            % Calculate Number of Padding Zeroes
            if any(strcmp(targetFilterType,{'bandpass','highpass'}))
                bpo = ceil(1/log2(f2/f1));
                nPeriods = 11.5*bpo;
                Lz = nPeriods * fs/f2 + 1230; % length of zero-padded signal (single-sided padding)    
            else % 'lowpass'
                nPeriods = 30; % TEST!!
                Lz = nPeriods * fs/f2 + 1230; % length of zero-padded signal (single-sided padding)
            end
            Lz = 2^ceil(log2(Lz)); % round to highest power of 2
            Lz = max(Lz,L);
            nZeros = Lz - L; % number of zeros for single-sided zero padding

            % Apply Zero Padding
            if strcmp(filtMode,'filter')
                x = [x; zeros(nZeros,1)];
            else % filtfilt
                Lz = L + 2*nZeros; % length of zero-padded signal (double-sided padding)
                x = [zeros(nZeros,1); x; zeros(nZeros,1)]; % double-sided zero-padding
            end
        end

       % Decimate
        nDecimSteps = log2(decimFactor);
        for n = 1:nDecimSteps % decimate X by 2 nDecimSteps times
            x = myfiltfilt(decimFilter,x);
            x = downsample(x,2);
        end

        % Filter
        switch metricsOutput
            case 0 % Calculate Filtered Waveform at Original Sampling Rate
                xf = filtFun(targetFilter,x);
                xbi = xf;
                for p = 1:nDecimSteps
                    Li = ceil(Lz/2^(nDecimSteps-p)); % target size of current xbi vector
                    xbi = 2 * upsample(xbi,2); % upsample by 2
                    % Trim or Pad signal to Right Length
                    if strcmp(filtMode,'filter')
                        nZeros = max(Li - length(xbi),0);
                        nTrimSamples = max(length(xbi) - Li,0);
                        startSample = 1;
                        endSample = length(xbi) - nTrimSamples;
                        xbi = [xbi(startSample:endSample); ...
                            zeros(nZeros,1)]; % interpolated signal trimmed to right size
                    else % filtfilt
                        nZeros = max(Li - length(xbi),0);
                        nZerosLeft = round(nZeros/2);
                        nZerosRight = nZeros - nZerosLeft;
                        nTrimSamples = length(xbi) - Li;
                        nTrimSamplesLeft = round(nTrimSamples/2);
                        nTrimSamplesRight = nTrimSamples - nTrimSamplesLeft;
                        startSample = 1 + nTrimSamplesLeft;
                        endSample = length(xbi) - nTrimSamplesRight;
                        xbi = [zeros(nZerosLeft,1); ...
                            xbi(startSample:endSample); ...
                            zeros(nZerosRight,1)]; % interpolated signal trimmed/padded to right size    
                    end
                    % Decimate
                    xbi = filtFun(decimFilter,xbi); % lowpass-filter upsampled signal
                end
                % Data Wrapping
                if dataWrap && zeroPad
                    xbw = datawrap(xbi,L); % warp over original length
                    if strcmp(filtMode,'filtfilt')
                        nZeros = round((Lz - L)/2); % zero-padding samples (one-sided) for original signal
                        Lc = rem(nZeros,L); % number of samples for circular shift
                        xbw = circshift(xbw,-Lc);
                    end
                    weight = sqrt(sum(xbi.^2)/sum(xbw.^2)); % correction factor to compensate for any energy cancellation from wrapping
                    xbi = xbw*weight;
                end
                xb = xbi; % filtered signal
            case 1 % Calculate Band Metrics of Input Signal x
                weight = sqrt(Lz/L); % correction factor to compensate the effect of zero-padding on the RMS
                xf = filtFun(targetFilter,x); % filter signal with target filter
                xb(1) = rms(xf)*weight; % root-mean-square value
                xb(2) = xb(1)^2*L/fs; % exposure value
                xb(3) = max(abs(xf)); % peak value
                xb(4) = max(xf) - min(xf); % peak-to-peak value
        end
    end
else
    switch metricsOutput
        case 0
            xb = x;
        case 1
            xb(1) = rms(x); % root-mean-square value
            xb(2) = xb(1)^2*L/fs; % exposure value
            xb(3) = max(abs(x)); % peak value
            xb(4) = max(x) - min(x); % peak-to-peak value
    end
end
