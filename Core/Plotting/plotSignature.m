% h = plotSignature(Radiation,varargin)
%
% DESCRIPTION
% Acoustic signature of the array calculated at the specified azimuth angle
% ('TargetAzimuth') and elevation angle ('TargetElevation'). The signature is 
% extracted from the RADIATION structure. The position of the primary and 
% bubble peaks is shown with circles. The signature can reflect the actual
% amplitude at the receiver location or be normalised to 1 m based on input 
% property 'Backprop'.
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
% - 'Backprop': TRUE if the amplitude of the signature is to be normalised 
%    to 1 m and FALSE otherwise (TRUE by default).
%
% OUTPUT ARGUMENTS
% - h: handle to array signature figure.
%
% FUNCTION CALL
% 1. h = plotSignature(Radiation)
% 2. h = plotSignature(...,PROPERTY,VALUES)
%    Properties: 'targetAzimuth', 'targetElevation', 'Backprop'
% 
% FUNCTION DEPENDENCIES
% - plotSignature_ec
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% See also PLOTHORIZONTALDIRECTIVITY, PLOTVERTICALDIRECTIVITY, PLOTSPECTRUM

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function h = plotSignature(Radiation,varargin)

% ERROR CONTROL
[targetAzimuth,targetElevation,backprop] = ...
    plotSignature_ec(Radiation,varargin);

% FIND POINT
sigma_all = [Radiation.Signature.souToRecAzimuth]';
theta_all = [Radiation.Signature.souToRecElevation]';
[~,iSigma] = min(abs(targetAzimuth - sigma_all));
[~,iTheta] = min(abs(targetElevation - theta_all));
targetAzimuth = sigma_all(iSigma); % exact target azimuth
targetElevation = theta_all(iTheta); % exact target elevation
iPoint = find(sigma_all == targetAzimuth & theta_all == targetElevation);

% LOAD PARAMETERS
fs = Radiation.Config.sampleRate;
pbr = Radiation.Metrics(iPoint).primaryToBubble; % primary to bubble ratio
bp = Radiation.Metrics(iPoint).bubblePeriod *1e3; % bubble period [ms]
primaryPeak = Radiation.Metrics(iPoint).primaryPeakPositive;
primaryTime = Radiation.Metrics(iPoint).primaryPeakPosTime *1e3; % main peak time [ms]
bubblePeak = Radiation.Metrics(iPoint).bubblePeakPositive;
bubbleTime = Radiation.Metrics(iPoint).bubblePeakPosTime *1e3; % bubble peak time [ms]

% LOAD SIGNATURE
r = Radiation.Signature(iPoint).souToRecDistance;
p = Radiation.Signature(iPoint).pressureAtOneMetre;
p = p(:); % convert pressure into column vector
nSamples = length(p);
t = (0:nSamples-1)/fs * 1e3; % time axis [ms]

% UPDATE SPECTRUM ('backprop')
if ~backprop
    p = p/r;
end

% UNITS
unitStr = 'Pa';
if backprop
    unitStr = strcat(unitStr,' @1m');
else
    unitStr = strcat(unitStr,sprintf(' @%0.0fm',r));  
end
    
% PLOT SIGNATURE
h = figure;
hold on
hplo(1) = plot(t,p,'Color','b','LineWidth',1.5);
yLim = get(gca,'YLim');
fill([-50 -50 0 0],[yLim(1) yLim(2) yLim(2) yLim(1)],[0.97 0.97 0.97],...
    'EdgeColor','none')
hplo(2) = plot(primaryTime,primaryPeak,'ro','MarkerSize',5);
hplo(3) = plot(bubbleTime,bubblePeak,'co','MarkerSize',5);
xlabel('Time [ms]')
ylabel(sprintf('Amplitude [%s]',unitStr))
title({'Acoustic Pressure Signature of Air Gun Array'; ...
    sprintf(['\\rm\\sigma = %0.1f°, \\theta = %0.1f°, PBR = %0.1f, '...
    '\\tau_b = %0.0f ms'],targetAzimuth,targetElevation,pbr,bp);''})
legend(hplo,{'Signature','Primary','Bubble'},'Location','NorthEast')
xlim([-50 max(t)])
ylim(yLim)
pbaspect([1 1 1])
set(gcf,'units','normalized','outerposition',[0.3 0.1 0.45 0.8])
box on
