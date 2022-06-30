% Metrics = signatureMetrics(x,fs,freqLimits)
%
% DESCRIPTION
% Computes several metrics from the input acoustic signature X with a sampling
% rate FS. The metrics include the band and broadband RMS/SEL amplitudes
% calculated within the frequency range given by FREQLIMITS, the time and
% magnitude of the positive and negative peaks of the bubble and primary pulse,
% the peak-to-peak amplitude of the bubble and primary pulse, the peak to
% bubble ratio (PBR) and the bubble period.
%
% SIGNATUREMETRICS.m calls FFTBANKFILTERDESIGN.m and FFTBANKFILTER.m to design
% a FFT filterbank and obtain the RMS/SEL metrics for the third-octave bands
% within the specified frequency range.
%
% INPUT ARGUMENTS
% - x: acoustic signature
% - fs: sampling rate [Hz]
% - freqLimits: two-element numeric vector containing the frequency limits
%   [fmin fmax] to use for processing the metrics.
%
% OUTPUT ARGUMENTS
% - Metrics: structure containing the metrics related to the input acoustic 
%   signature X. For details about the content of this structure see funtion
%   INITIALISEMETRICSSTRUCT.m
%
% FUNCTION CALL
% 1. Metrics = signatureMetrics(x,fs,freqLimits)
% 
% FUNCTION DEPENDENCIES
% - fftBankFilterDesign
% - fftBankFilter
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% See also BUBBLERADIATION

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function Metrics = signatureMetrics(x,fs,freqLimits)

% Filter Design
bpo = 3; % bands per octave (bpo = 3 for third-octave)
FftFilterBank = fftBankFilterDesign(fs,bpo,freqLimits);
fc = FftFilterBank.centralFreq; % central frequencies [Hz]
fn = FftFilterBank.nominalFreq; % nominal central frequencies [Hz]
zeroPad = true;

% Main Peak
[z2p_pri,iz2p_pri] = max(x); % zero-to-positive peak (amplitude and index)
[z2n_pri,iz2n_pri] = min(x); % zero-to-negative peak (amplitude and index)
p2p_pri = z2p_pri - z2n_pri; % peak-to-peak (amplitude)
tz2p_pri = (iz2p_pri - 1)/fs;
tz2n_pri = (iz2n_pri - 1)/fs;

% Bubble
iZeroCross = find(x(iz2p_pri:end) < 0,1,'first') + iz2p_pri - 1;
[~,ind] = max(x(iZeroCross:end));
iz2p_bub = ind + iZeroCross - 1;
[~,ind] = min(x(iz2p_bub :end));
iz2n_bub = ind + iz2p_bub - 1;
z2p_bub = x(iz2p_bub);
z2n_bub = x(iz2n_bub);
p2p_bub = z2p_bub - z2n_bub; 
tz2p_bub = (iz2p_bub - 1)/fs;
tz2n_bub = (iz2n_bub - 1)/fs;

% Signature Parameters
pbr = z2p_pri/z2p_bub;
T_bub = tz2p_bub - tz2p_pri;

% Root Mean Square
xb = fftBankFilter(FftFilterBank,x,zeroPad);
x_rmsb = xb(1,:);
x_selb = xb(2,:);
x_rms = sqrt(sum(x_rmsb.^2));
x_sel = sum(x_selb);

% Populate Metrics Structure
Metrics.primaryToBubble = pbr;
Metrics.bubblePeriod = T_bub;
Metrics.primaryPeakPositive = z2p_pri;
Metrics.primaryPeakNegative = z2n_pri;
Metrics.primaryPeakToPeak = p2p_pri;
Metrics.primaryPeakPosTime = tz2p_pri;
Metrics.primaryPeakNegTime = tz2n_pri;
Metrics.bubblePeakPositive = z2p_bub;
Metrics.bubblePeakNegative = z2n_bub;
Metrics.bubblePeakToPeak = p2p_bub;
Metrics.bubblePeakPosTime = tz2p_bub;
Metrics.bubblePeakNegTime = tz2n_bub;
Metrics.signatureRmsBroadband = x_rms;
Metrics.signatureSelBroadband = x_sel;
Metrics.signatureRmsBand = x_rmsb;
Metrics.signatureSelBand = x_selb;
Metrics.centralFrequencies = fc;
Metrics.nominalFrequencies = fn;

