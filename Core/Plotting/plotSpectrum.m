% h = plotSpectrum(Radiation,varargin)
%
% DESCRIPTION
% Frequency spectrum of the acoustic signature of the array calculated at the
% specified azimuth ('TargetAzimuth') and elevation ('TargetElevation') angles.
% The signature is extracted from the RADIATION structure. The Power Spectral
% Density ('psd') or the band spectrum ('cpb') can be represented based on
% the input property 'spectrumType'. The property 'Metric' determines whether
% the root-mean-square ('rms') or the sound exposure ('sel') levels are used
% in the spectral representation. The spectrum can reflect the actual amplitude
% at the receiver location or be normalised to 1 m based on input property 
% 'Backprop'.
%
% INPUT ARGUMENTS
% - Radiation: structure containing the acoustic signature of the array at
%   one or more receiver locations, related metrics, and configuration 
%   information. For details about the content of this structure see funtion
%   INITIALISERADIATIONSTRUCT.m
%
% INPUT PROPERTIES
% - 'TargetAzimuth': azimuth between the centre of the array and the receiver
%   [deg]. The angle is given in geographic format. The value must be between
%   -180 and 180 deg, where 90 represents the North direction. A value of 
%   TARGETAZIMUTH = 0 deg is used by default.
% - 'TargetAzimuth': elevation between the centre of the array and the receiver
%   [deg]. The angle is given in geographic format. The value must be between
%   0 and 90 deg, where 90 represents the vertical direction. A value of 
%   TARGETELEVATION = 90 deg is used by default.
% - 'SpectrumType': type of spectral representation.
%   # 'psd': power spectral density spectrum
%   # 'cpb': constant percentual band spectrum
% - 'Metric': acoustic metric for the spectrum
%   # 'rms': root-mean-square
%   # 'sel': sound exposure level
% - 'Backprop': TRUE if the amplitude of the signature is to be normalised 
%    to 1 m and FALSE otherwise (TRUE by default).
%
% OUTPUT ARGUMENTS
% - h: handle to the frequency spectrum figure.
%
% FUNCTION CALL
% 1. h = plotSpectrum(Radiation)
% 2. h = plotSpectrum(...,PROPERTY,VALUES)
%    Properties: 'targetAzimuth', 'targetElevation', 'SpectrumType',
%       'Metric','Backprop'
% 
% FUNCTION DEPENDENCIES
% - plotSpectrum_ec
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% See also PLOTHORIZONTALDIRECTIVITY, PLOTVERTICALDIRECTIVITY, PLOTSPECTRUM

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function h = plotSpectrum(Radiation,varargin)

% ERROR CONTROL
[targetAzimuth,targetElevation,spectrumType,metric,backprop] = ...
    plotSpectrum_ec(Radiation,varargin);

% LOAD PARAMETERS
fs = Radiation.Config.sampleRate;
pref = 1e-6; % reference pressure [uPa]

% FIND POINT
sigma_all = [Radiation.Signature.souToRecAzimuth]';
theta_all = [Radiation.Signature.souToRecElevation]';
[~,iSigma] = min(abs(targetAzimuth - sigma_all));
[~,iTheta] = min(abs(targetElevation - theta_all));
targetAzimuth = sigma_all(iSigma); % exact target azimuth
targetElevation = theta_all(iTheta); % exact target elevation
iPoint = find(sigma_all == targetAzimuth & theta_all == targetElevation);

% LOAD SIGNATURE
r = Radiation.Signature(iPoint).souToRecDistance;
p = Radiation.Signature(iPoint).pressureAtOneMetre;
p = p(:); % convert pressure into column vector
nSamples = length(p);

% COMPUTE RMS SPECTRUM
if strcmpi(spectrumType,'psd')
    nFft = 2^14;
    P = fft(p,nFft);
    P = P .* conj(P)/(nFft*fs);
    P = P * nFft/nSamples; % correct effect of zero-padding
    iLast = round(nFft/2);
    P = [P(1); 2*P(2:iLast)];
    f = (0:iLast-1)*fs/(2*iLast);
    
else % spectrumType = 'cpb'
    f = Radiation.Metrics(iPoint).centralFrequencies;
    fn = Radiation.Metrics(iPoint).nominalFrequencies;
    P = Radiation.Metrics(iPoint).signatureRmsBand.^2;  
end

% UPDATE SPECTRUM ('metric')
if strcmpi(metric,'sel')
    P = P * nSamples/fs;
end

% UPDATE SPECTRUM ('backprop')
if ~backprop
    P = P/r;
end

% UNITS
if strcmpi(spectrumType,'psd')
    if strcmpi(metric,'rms')
        unitStr = 'Pa^2/Hz';
    else % metric = 'sel'
        unitStr = 'Pa^2s/Hz';
    end
else % spectrumType = 'cpb'
    if strcmpi(metric,'rms')
        unitStr = 'Pa^2';
    else % metric = 'sel'
        unitStr = 'Pa^2s';
    end 
end
if backprop
    unitStr = strcat(unitStr,' @1m');
else
    unitStr = strcat(unitStr,sprintf(' @%0.0fm',r));  
end
    
% PLOT SPECTRUM
if strcmpi(spectrumType,'psd')
    h = figure;
    P_db = 10*log10(P/pref^2);
    plot(f,P_db,'Color','b','LineWidth',1.5)
    xlabel('Frequency [Hz]')
    ylabel(sprintf('PSD [dB re %s]',unitStr))
    title({'Power Spectral Density Spectrum'; sprintf(['\\rm\\sigma '...
        '= %0.1f°, \\theta = %0.1f°, %s'],targetAzimuth,targetElevation,...
        upper(metric));''})
    set(gca,'XScale','log')
    ymax = ceil(max(P_db)/10)*10;
    ymin = ymax - 100;
    ylim([ymin ymax])
    xlim([10 max(f)])
    grid on
    pbaspect([1 1 1])
    set(gcf,'units','normalized','outerposition',[0.3 0.1 0.45 0.8])
else % spectrumType = 'cpb'
    h = figure;
    P_db = 10*log10(P/pref^2);
    plot(f,P_db,'Color','b','Marker','o','MarkerFaceColor','b',...
        'MarkerSize',3,'LineWidth',1.5)
    xlabel('Frequency [Hz]')
    ylabel(sprintf('CPB [dB re %s]',unitStr))
    title({'Third-Octave Band Spectrum'; sprintf(['\\rm\\sigma = '...
        '%0.1f°, \\theta = %0.1f°, %s'],targetAzimuth,targetElevation,...
        upper(metric));''})
    set(gca,'XScale','log')
    iMajorTicks = false(length(fn),1);
    iMajorTicks(1:3:length(fn)) = true;
    xMajorTicks = fn(iMajorTicks);
    xMinorTicks = fn(~iMajorTicks);
    set(get(gca,'XAxis'),'MinorTickValues',xMinorTicks)
    f_name = compactnum(fn(1:3:end),1);
    set(gca,'XTick',xMajorTicks)
    set(gca,'XTickLabel',f_name)
    ymax = ceil(max(P_db)/10)*10;
    ymin = ymax - 50;
    ylim([ymin ymax])
    xlim([min(f) max(f)])
    grid on
    pbaspect([1 1 1])
    set(gcf,'units','normalized','outerposition',[0.3 0.1 0.45 0.8])
end

