%  [f,fcnom] = FRACTIONALOCTAVEBANDS(varargin)
%
%  DESCRIPTION
%  Calculates the central and bandedge frequencies of the fractional octave 
%  bands that are within the selected frequency range. The calculations are 
%  based on the BS EN 61260-1:2014 standard (equivalent to ANSI S1.11:2014).
%
%  INPUT VARIABLES (Properties)
%  - 'BandsPerOctave': bandwidth factor or number of bands in an octave 
%    (i.e. BPO = 1, octave; BPO = 2, half-octave; BPO = 3, third-octave, etc). 
%    Reciprocal of the bandwidth designator (see "Terms and Definitions" 
%    in BS EN 61260-1:2014 for definition).
%  - 'FrequencyLimits': 2 element vector containing the frequency limit
%    for the calculations [fmin fmax]. The function returns the exact and
%    nominal frequencies of those band-pass filters whose side frequencies 
%    fall within the specified range.
%  - 'LimitMode': string describing the singular point in the bandpass
%    filter at which the limit frequencies are defined. There options:
%    ¬ 'central': the frequency limits refer to the central point of the 
%      top and bottom bands (fmin = fcmin, fmax = fcmax)
%    ¬ 'bandedge': the frequency limits refer to the higher and lower -6dB 
%      point of the top and bottom bands, respectively (fmin = fcmin *
%      G^(-0.5/BPO), fmax = fcmax * G^(0.5/BPO)).
%    ¬ 'bandstop': the frequency limits refer to the higher and lower
%      band stop point of the top and bottom bands, respectively (fmin = 
%      fcmin * G^(-1/BPO), fmax = fcmax * G(1/BPO)).
%  - 'Base10': flag specifying the base G of the band factor G^(1/BPO).
%    ¬ 0: base 2 band factor, i.e. G = 2 (e.g. 1000 Hz reference, 1000*2 
%      higher octave). Avoid
%    ¬ 1: base 10 band factor, i.e. G = 10^(3/10) (e.g. 1000 Hz reference,
%      1000*10^(3/10) higher octave). This is the DEFAULT option.
%
%  OUTPUT VARIABLES
%  - f: three column array including the exact frequencies for the lower
%    band edge f1, band centre fc and higher band edge f2 ([f1 fc f2]) [Hz]
%  - fnom: column vector including the nominal central frequencies [Hz]
%
%  COMMENTS
%  - A base-2 calculation should be avoided, as this will result in frequency
%    inaccuracies. The error will be larger that further away the frequency
%    is from the reference of 1 kHz (see BS EN 61260-1:2014, Annex D). These
%    inaccuracies also affect the calculation of nominal frequencies for
%    B >> 24.
%
%  - The Rounding Method on Narrow-Band Filters: as the filter bandwidth
%    becomes narrower (i.e. higher BPO number) the central frequencies
%    become more similar. If the rounding method is applied to the same
%    number of significant digits used for 1-octave band filters it may
%    occur that contiguous filters share the same nominal frequency. To 
%    avoid this situation, the number of digits used for rounding should be
%    increased with B. 
%
%  - A reasonable approach would be to increase by one the number of rounding 
%    digits when two contiguous central frequencies differ less than 10%, 
%    increase by two rounding digits when the difference is less than 1%, 
%    and so on. The ratio of central frequencies from consecutive bands is 
%    given by G^(1/B), where G is the one-octave factor (~2); when G^(1/B) - 1
%    < 0.1 the number of rounding digits should be increased by 1. The 
%    following table shows the increase in rounding digits with the nominal 
%    filter bandwidth.
%   
%  - In this function we stick to the standard by not increasing the
%    number of significant digits to be used for rounding. Therefore, 
%    calculations of nominal frequencies will not be accurate for B >> 24.
%
%    INCREASE IN ROUNDING DIGITS Ddig WITH FILTER BANDWIDTH
%
%       | BPO | G^(1/BPO) | Ddig |
%     ----------------------------
%          1      2       0
%          2      1.41    0
%          3      1.26    0
%        ...      ...     ...
%          6      1.12    0
%          7      1.104   0
%          8      1.09    1
%        ...      ...     ...
%         24
%         70      1.0099  2
%     ----------------------------    
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  FUNCTION CALLS
%  - [f,fcnom] = fractionalOctaveBands(property,value)  
%    ¬ properties: 'BandsPerOctave', 'FrequencyLimits','LimitMode','Base10'

%  VERSION 2.0
%  Revised by: Guillermo Jimenez
%  Date: 20 May 2021
%  - Modified calculation method for nominal central frequencies to stick
%    to the BS EN standard.
%
%  VERSION 1.1 
%  Revised by: Guillermo Jimenez
%  Date: 20 May 2019
%  - Call modified to allow selection of parameters by property/value pair.
%  - New property 'LimitMode' to define the singular point that the limit
%    frequencies refer to ('central','bandedge','bandstop'). Before, the
%    option equivalent to 'bandedge' was the only available.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  17 Jan 2018

function [f,fcnom] = fractionalOctaveBands(varargin)

% Check Number of Input Arguments
M = nargin;
narginchk(0,8) % error if number of input arguments is higher than 6
if rem(M,2)
    error('Property and value input arguments must come in pairs')
end

% Initialise Default Parameters
bpo = 1;
frange = [20 20e3];
fmode = 'central';
base10 = 1;

% Retrieve Input Variables
for m = 1:2:M
    inputProperty = varargin{m};
    if ~ismember(inputProperty,{'BandsPerOctave','FrequencyLimits','LimitMode','Base10'})
        error('Invalid input property')
    else
        switch inputProperty
            case 'BandsPerOctave'
                bpo = varargin{m+1};
            case 'FrequencyLimits'
                frange = varargin{m+1};
            case 'LimitMode'
                fmode = varargin{m+1};
            case 'Base10'
                base10 = varargin{m+1};
        end
    end
end

% Set Constant based on Selected LimitMode
fmode = lower(fmode);
switch fmode
    case 'central'
        k = 0;
    case 'bandedge'
        k = 1/2;
    case 'bandstop'
        k = 1;    
    otherwise
        error('Value for LIMITMODE input property not recognised')
end
    
% Calculation of Exact Frequencies ¦f1¦, ¦fc¦, ¦f2¦
if base10, G = 10^(3/10); else, G = 2; end % one octave frequency ratio
fref = 1000; % reference frequency of 1 kHz
fmin = frange(1); % lower frequency limit [Hz]
fmax = frange(2); % higher frequency limit [Hz]
a = double(~rem(bpo,2))/2; % 0 for BPO = odd, 1/2 for BPO = even
xmax = floor(bpo*log(fmax/fref)/log(G) - k - a); % maximum band index (xmax = no. bands above the reference band of 1 kHz)
xmin = ceil(bpo*log(fmin/fref)/log(G) + k - a); % minimum band index (abs(xmin) = no. of bands below the reference band of 1 kHz)
x = xmin:xmax; % vector of band indices
fc = fref * G.^((x + a)/bpo)'; % vector of exact central frequencies from selected bands [Hz]
f1 = fc * G^(-1/(2*bpo)); % vector of exact lower frequencies from selected bands [Hz]
f2 = fc * G^(1/(2*bpo)); % vector of exact higher frequencies from selected bands [Hz] 
f = [f1 fc f2]; % array of exact frequencies

% Calculation of Nominal Central Frequencies ¦fcnom¦
if any(bpo == [1 3])
    fcn = [1, 1.25, 1.6, 2, 2.5, 3.15, 4, 5, 6.3, 8, 10]; % third-octave central frequencies normalised to 1 MSD
    ndig = floor(log10(fc)) + 1; % number of digits in the integer part of the decimal central frequencies
    [~,idx] = min(abs(fc - fcn.*10.^(ndig-1)),[],2);
    fcnom = fcn(idx)'.*10.^(ndig-1);    
elseif bpo == 2
    fcnom = round(fc,3,'significant');
else % bpo > 3
    msd = floor(fc.*10.^(-floor(log10(fc)))); % vector with the Most Significant Digit (MSD) of each central frequency
    idx_3dig = msd < 5; % logical vector indicating the position in ¦fc¦ of those frequencies with the MSD = 1 to 4
    idx_2dig = msd >= 5; % logical vector indicating the position in ¦fc¦ of those frequencies with the MSD = 5 to 9
    fcnom_3dig = round(fc,3,'significant');
    fcnom_2dig = round(fc,2,'significant');
    fcnom(idx_3dig,1) = fcnom_3dig(idx_3dig); 
    fcnom(idx_2dig,1) = fcnom_2dig(idx_2dig);  
end