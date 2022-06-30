% h = plotVerticalDirectivity(Radiation,varargin)
%
% DESCRIPTION
% Vertical directivity of the array calculated at the specified azimuth
% angle ('TargetAzimuth') and frequencies ('Frequencies') based on the
% acoustic signature data contained in the RADIATION structure.
%
% INPUT ARGUMENTS
% - Radiation: structure containing the acoustic signature of the array at
%   one or more receiver locations, related metrics, and configuration 
%   information. For details about the content of this structure see funtion
%   INITIALISERADIATIONSTRUCT.m
%
% INPUT PROPERTIES
% - 'TargetAzimuth': azimuth at which the vertical plane crosses the 3D 
%    emission pattern of the air gun array [deg]. The angle is given in
%    geographic format. The value must be between -180 and 180 deg, with 90
%    representing the North direction. A value of TARGETAZIMUTH = 0 deg
%    is used by default.
% - 'Frequencies': vector of frequencies at which the directivity pattern
%    is computed [Hz]. The frequencies in RADIATION closest to the specified
%    will be used. For broadband directivity use FREQUENCIES = [] (default)
%
% OUTPUT ARGUMENTS
% - h: handle to horizontal directivity figure.
%
% FUNCTION CALL
% 1. h = plotVerticalDirectivity(Radiation)
% 2. h = plotVarticalDirectivity(...,PROPERTY,VALUE)
%    Properties: 'targetAzimuth', 'Frequencies'
% 
% FUNCTION DEPENDENCIES
% - plotVerticalDirectivity_ec
% - below360
% - semiPolarPlot
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% See also PLOTHORIZONTALDIRECTIVITY, PLOTSIGNATURE, PLOTSPECTRUM

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function h = plotVerticalDirectivity(Radiation,varargin)

% ERROR CONTROL
[targetAzimuth,frequencies] = ...
    plotVerticalDirectivity_ec(Radiation,varargin);

% FIND FREQUENCIES
fc_orig = [Radiation.Metrics(1).centralFrequencies];
fc = unique(interp1(fc_orig,fc_orig,frequencies,'nearest'));
iFreqs = find(ismember(fc_orig,fc));
nFreqs = length(iFreqs);

% FIND POINTS
% All Angles
sigma_all = below360([Radiation.Signature.souToRecAzimuth],'deg+');
theta_all = below360([Radiation.Signature.souToRecElevation],'deg+');

% Points in Right Quadrant
targetAzimuth1 = below360(targetAzimuth,'deg+');
[~,ind] = min(abs(sigma_all - targetAzimuth1));
targetAzimuth1 = sigma_all(ind); % exact target elevation in data
i1 = find(sigma_all == targetAzimuth1);

% Vertical Point
targetElevation = 90;
[~,ind] = min(abs(theta_all - targetElevation));
targetElevation = theta_all(ind);
i2 = find(theta_all == targetElevation);

% Points in Left Quadrant
targetAzimuth2 = below360(targetAzimuth + 180,'deg+');
[~,ind] = min(abs(sigma_all - targetAzimuth2));
targetAzimuth2 = sigma_all(ind); % exact target elevation in data
i3 = find(sigma_all == targetAzimuth2);

% PROCESS DIRECTIVITY
if isempty(frequencies) % broadband
    % Elevation Angles (math form: 0 E, 90 N)
    theta1 = -[Radiation.Signature(i1).souToRecElevation]'; % right quadrant
    theta2 = -[Radiation.Signature(i2).souToRecElevation]'; % centre point 
    theta3 = [Radiation.Signature(i3).souToRecElevation]' - 180; % left quadrant
    
    % Directivity
    directivity1 = [Radiation.Metrics(i1).signatureRmsBroadband]'; % right quadrant
    directivity2 = [Radiation.Metrics(i2).signatureRmsBroadband]'; % centre point
    directivity3 = [Radiation.Metrics(i3).signatureRmsBroadband]'; % left quadrant
    
    % Sort Parts
    [~,iSort1] = sort(theta1,'descend');
    [~,iSort3] = sort(theta3,'descend');
    theta1 = theta1(iSort1);
    theta3 = theta3(iSort3);
    directivity1 = directivity1(iSort1,:);
    directivity3 = directivity3(iSort3,:);
    
    % Final Angles and Directivity
    theta = [theta1; theta2; theta3];
    directivity = [directivity1; directivity2; directivity3];
    directivity_norm = directivity/max(directivity);
    
    % Legend
    legendText = {'Broadband'};
    
else % band
    % Elevation Angles (math form: 0 E, 90 N)
    theta1 = -[Radiation.Signature(i1).souToRecElevation]'; % right quadrant
    theta2 = -[Radiation.Signature(i2).souToRecElevation]'; % centre point
    theta3 = [Radiation.Signature(i3).souToRecElevation]' - 180; % right quadrant
    
    % Directivity
    directivity1 = reshape([Radiation.Metrics(i1).signatureRmsBand],...
        length(fc_orig),length(theta1))'; % right quadrant
    directivity2 = [Radiation.Metrics(i2).signatureRmsBand]; % centre point
    directivity3 = reshape([Radiation.Metrics(i3).signatureRmsBand],...
        length(fc_orig),length(theta3))'; % left quadrant
    
    % Sort Parts
    [~,iSort1] = sort(theta1,'descend');
    [~,iSort3] = sort(theta3,'descend');
    theta1 = theta1(iSort1);
    theta3 = theta3(iSort3);
    directivity1 = directivity1(iSort1,:);
    directivity3 = directivity3(iSort3,:);
    
    % Final Angles and Directivity
    theta = [theta1; theta2; theta3];
    directivity = [directivity1(:,iFreqs); directivity2(iFreqs); ...
        directivity3(:,iFreqs)];
    directivity_norm = directivity./max(directivity);
    
    % Legend
    for m = 1:nFreqs
        legendText{m} = sprintf('%0.0f Hz',fc(m));
    end
end

% Interpolate Directivity
nInterp = 180; % number of interpolation points
interpFactor = max(round(nInterp/length(theta)),1);
interpPoints = (1:1/interpFactor:length(theta))';
theta_int = interp1(theta,interpPoints,'linear');
directivity_int = interp1(directivity_norm,interpPoints,'pchip');
    
% PLOT VERTICAL DIRECTIVITY
theta_plo = [theta_int; theta_int(1)];
directivity_plo = [directivity_int; directivity_int(1,:)];
h = semiPolarPlot(theta_plo,directivity_plo,[0 1],legendText);
title({'Vertical Directivity of Air Gun Array';...
    sprintf('\\rm\\sigma = %0.0f to %0.0f (L-R)',...
    below360(targetAzimuth+180,'deg'),below360(targetAzimuth,'deg'));''})

