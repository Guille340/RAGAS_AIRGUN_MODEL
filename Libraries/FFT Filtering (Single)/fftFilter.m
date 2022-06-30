%  xb = FFTFILTER(x,fn1,fn2,varargin)
%
%  DESCRIPTION
%  Calculates the RMS amplitude of X within the frequency band delimited by 
%  the normalised cutoff frequencies FN1 and FN2. 
%
%  The variable input argument ZEROPAD is a logical value that indicates 
%  whether zero padding will be applied before filtering. Applying zero-
%  padding is particularly important when the period associated to the 
%  frequency band to be filtered is comparable to or larger than the duration 
%  of the filtered signal X. ZEROPAD = TRUE is used by default. As a general
%  rule, signals should not be filtered at frequencies lower than 10/TAU,
%  where TAU is the duration of the signal.
%
%  INPUT ARGUMENTS
%  - x: signal to be filtered.
%  - fn1: normalised lower cutoff frequency of the bandpass filter.
%  - fn2: normalised upper cutoff frequency of the bandpass filter.
%  - zeroPad (varargin{1}): logical value indicating whether zero-padding
%    will be applied. ZEROPAD = TRUE by default.
%
%  NOTE: FN1 and FN2 must be between 0 and 1, with 1 corresponding to
%  the Nyquist frequency (= sampleRate/2).
%
%  OUTPUT ARGUMENTS
%  - xb: vector of RMS values for each frequency band delimited by input
%    cutoff frequencies FN1 and FN2.
%
%  FUNCTION CALL
%  1. xb = FFTFILTER(x,fn1,fn2) -> zeroPad = TRUE
%  2. xb = FFTFILTER(x,fn1,fn2,zeroPad)
%
%  FUNCTION DEPENDENCIES
%  - None
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  CONSIDERATIONS AND LIMITATIONS
%  - For calculating the exposure amplitude XE use the formula below:
%    
%    XE = XB.^2 * LENGTH(X)/FS, where LENGTH(X)/FS is the signal's duration
%
%  See also FRACTIONALOCTAVEBANDS

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  26 Jul 2021

function xb = fftFilter(x,fn1,fn2,varargin)

narginchk(3,4)

% Default Variable Input Arguments
zeroPad = true;

% Assign Variable Input Arguments
if nargin == 4
    zeroPad = varargin{1};
end

% Error Control
if ~isnumeric(x)
    error('X must be a numeric vector or matrix')
end
if ~isnumeric(fn1) || ~isvector(fn1) || any(fn1 < 0 | fn1 > 1)
    error('FN1 must be a numeric vector with values between 0 and 1')
end
if ~isnumeric(fn2) || ~isvector(fn2) || any(fn2 < 0 | fn2 > 1)
    error('FN1 must be a numeric vector with values between 0 and 1')
end
if ~any(zeroPad == [0 1])
    zeroPad = true;
    warning(['ZEROPADDING property must be [0,1] or logical. '...
        'ZEROPADDING = TRUE will be used'])
end

% General Parameters
if isvector(x), x = x(:); end % make x a column vector
fnc = sqrt(fn1.*fn2); % normalised central frequencies
bpo = max(round(1./log2(fn2./fn1))); % bands per octave
nSamples = size(x,1); % original length of signals in X 
nSignals = size(x,2); % number of signals

% Zero-Padding    
nSamplesPad = nSamples; % initialise number of zero-padded samples
if zeroPad 
    nPeriods = 130*bpo^1.23;
    nSamplesPad = nPeriods * 1/(2*min(fnc));
    nSamplesPad = max(2^ceil(log2(nSamplesPad)),nSamples);
    nZeros = nSamplesPad - nSamples;
    x = [x; zeros(nZeros,nSignals)];
end

% Power Autospectrum (one-sided)
X = fft(x);
X = X.*conj(X)/nSamplesPad^2;
if rem(nSamplesPad,2) % odd number of samples
    X = [X(1,:); 2*X(2:(nSamplesPad+1)/2,:)]; 
else % even number of samples
    X = [X(1,:); 2*X(2:nSamplesPad/2,:); X(nSamplesPad/2 + 1,:)]; 
end
fn = (0:round(nSamplesPad/2))'*2/nSamplesPad; % normalised frequency vector (fn = 0:2/nSamplesPad:1)

% Calculation of Band RMS
nFraBands = length(fnc); % number of fractional octave bands
xb = nan(nFraBands,nSignals);
for k = 1:nFraBands
    isFreqsInBand = fn >= fn1(k) & fn < fn2(k);
    if ~isempty(isFreqsInBand)
        xb(k,:) = sqrt(sum(X(isFreqsInBand,:),1))*sqrt(nSamplesPad/nSamples); % rms values
    end
end