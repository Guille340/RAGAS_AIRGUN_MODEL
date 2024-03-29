%  xb = FFTSINGLEFILTER(FftFilter,x,varargin)
%
%  DESCRIPTION
%  Filters the input signal X with sampling rate FS using a FFT filter 
%  structure FFTFILTER. An appropriate amount of zero-padding is applied to
%  the signal by default (ZEROPAD = TRUE). The zero padding eliminates any 
%  issue that may arise when filtering signals with a length comparable to 
%  the filter�s group delay.
%
%  INPUT ARGUMENTS 
%  - FftFilter: structure containing information for the FFT digital.
%    FFTFILTER is generated with FFTSINGLEFILTERDESIGN (see function for
%    details about its fields).
%  - x: vector of data to be filtered. Its length must be three times larger, 
%    or more, than the filter order (if unknown, use filtord to extract the 
%    filter order).
%  - zeroPad (varargin{1}): logical or [0,1] numeric value indicating whether 
%    zero-padding will be applied to the input signal X. Zero-padding is
%    necessary in those cases where the signal is short enough that lowpass
%    or bandpass filtering will cause a delay that will displace the filtered 
%    signal beyond its original window length. This will result in a noticeable 
%    energy error. Depending on the signal duration and top cutoff frequency 
%    of the filter, the error can be as large as several tens of decibels. 
%    One-sided zero padding is applied when FILTERMODE = 'filter', and two-
%    sided zero padding when using FILTERMODE = 'filtfilt'. Any filtering 
%    approach will suffer from this "time leakage" effect. If this property 
%    is omitted, ZEROPADDING = TRUE (default).
%    
%  OUTPUT ARGUMENTS
%  - y: two-element vector containing the RMS and exposure amplitudes of
%    the filtered signal.
%
%  FUNCTION CALL
%  xb = FFTSINGLEFILTER(FftFilter,x) % ZEROPAD = TRUE
%  xb = FFTSINGLEFILTER(FftFilter,x,zeroPad)
%
%  FUNCTION DEPENDENCIES
%  - fftFilter
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  EXAMPLE
%  % 1) Configuration Data
%  fs = 44100;
%  freqCutoff = [20 10e3];
%  zeroPad = true;
%  T = 1; % signal duration [s]
%  x = 2*rand(1,T*fs) - 1; % signal
%
%  % 2) Design FFT Filter
%  FftFilter = fftSingleFilterDesign(fs,freqCutoff)
%
%  % 3) Filter Signal with FFT Method
%  xb_fft = fftSingleFilter(FftFilter,x,zeroPad)
%
%  % 4) Band Levels
%  Xb_fft_rms = 20*log10(xb_fft(1));
%  Xb_fft_exp = 10*log10(xb_fft(2));
%
%  See also FFTFILTER

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  07 Aug 2021

function xb = fftSingleFilter(FftFilter,x,varargin)

% Check Number of Input Arguments
narginchk(2,3)

% Assign Variable Input Arguments
zeroPad = true;
if nargin == 4
    zeroPad = varargin{1};
end

% Error Control (ARGUMENT: 'x')
if ~isnumeric(x) || ~isvector(x)
    error('Input argument X must be a numeric vector')
end

% Error Control (ARGUMENT: 'zeroPad')
if ~any(zeroPad == [0 1])
    zeroPad = true;
    warning(['ZEROPADDING property must be [0,1] or logical. '...
        'ZEROPADDING = TRUE will be used'])
end

% Load Filter Parameterws
fs = FftFilter.sampleRate;
fn1 = FftFilter.halfPowerFreqn1;
fn2 = FftFilter.halfPowerFreqn2;

% Filter Signal
xb_rms = fftFilter(x,fn1,fn2,zeroPad); % RMS of filtered signal
xb_exp = xb_rms.^2 * length(x)/fs; % exposure of filtered signal
xb(1) = xb_rms;
xb(2) = xb_exp;
