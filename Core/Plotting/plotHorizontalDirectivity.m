% h = plotHorizontalDirectivity(Radiation,varargin)
%
% DESCRIPTION
% Horizontal directivity of the array calculated at the specified elevation
% angle ('TargetElevation') and frequencies ('Frequencies') based on the
% acoustic signature data contained in the RADIATION structure.
%
% INPUT ARGUMENTS
% - Radiation: structure containing the acoustic signature of the array at
%   one or more receiver locations, related metrics, and configuration 
%   information. For details about the content of this structure see funtion
%   INITIALISERADIATIONSTRUCT.m
%
% INPUT PROPERTIES
% - 'TargetElevation': elevation at which the horizontal plane crosses the
%    3D emission pattern of the air gun array [deg]. The angle is given in
%    geographic format. The value must be between 0 and 90 deg, with 90
%    representing the vertical direction. A value of TARGETELEVATION = 45 deg
%    is used by default.
% - 'Frequencies': vector of frequencies at which the directivity pattern
%    is computed [Hz]. The frequencies in RADIATION closest to the specified
%    will be used. For broadband directivity use FREQUENCIES = [] (default)
%
% OUTPUT ARGUMENTS
% - h: handle to horizontal directivity figure.
%
% FUNCTION CALL
% 1. h = plotHorizontalDirectivity(Radiation)
% 2. h = plotHorizontalDirectivity(...,PROPERTY,VALUE)
%    Properties: 'targetElevation', 'Frequencies'
% 
% FUNCTION DEPENDENCIES
% - plotHorizontalDirectivity_ec
% - below360
% - polarPlot
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% See also PLOTVERTICALDIRECTIVITY, PLOTSIGNATURE, PLOTSPECTRUM

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function h = plotHorizontalDirectivity(Radiation,varargin)

% ERROR CONTROL
[targetElevation,frequencies] = ...
    plotHorizontalDirectivity_ec(Radiation,varargin);

h = []; % initialise figure handle
if targetElevation ~= 90
    % FIND FREQUENCIES
    fc_orig = [Radiation.Metrics(1).centralFrequencies];
    fc = unique(interp1(fc_orig,fc_orig,frequencies,'nearest'));
    iFreqs = find(ismember(fc_orig,fc));
    nFreqs = length(iFreqs);

    % FIND POINTS
    % All Points
    theta_all = below360([Radiation.Signature.souToRecElevation],'deg+');

    % Selected Points
    targetElevation = below360(targetElevation,'deg+');
    [~,ind] = min(abs(theta_all - targetElevation));
    targetElevation = theta_all(ind); % exact target elevation in data
    i1 = find(theta_all == targetElevation);

    % PROCESS DIRECTIVITY
    if isempty(frequencies) % broadband
        % Azimuth Angles (math form: 0 E, 90 N)
        sigma = 90 - [Radiation.Signature(i1).souToRecAzimuth]';

        % Directivity
        directivity = [Radiation.Metrics(i1).signatureRmsBroadband]';

        % Sort Values
        [~,iSort] = sort(sigma,'ascend');
        sigma = sigma(iSort);
        directivity = directivity(iSort);

        % Normalised Directivity
        directivity_norm = directivity/max(directivity);

        % Legend
        legendText = {'Broadband'};

    else % band
        % Azimuth Angles (math form: 0 E, 90 N)
        sigma = 90 - [Radiation.Signature(i1).souToRecAzimuth]';

        % Directivity
        directivity = reshape([Radiation.Metrics(i1).signatureRmsBand],...
            length(fc_orig),length(sigma))';

        % Sort Values
        [~,iSort] = sort(sigma,'ascend');
        sigma = sigma(iSort);
        directivity = directivity(iSort,:);

        % Normalised Directivity for Selected Frequencies
        directivity = directivity(:,iFreqs);
        directivity_norm = directivity./max(directivity);

        % Legend
        for m = 1:nFreqs
            legendText{m} = sprintf('%0.0f Hz',fc(m));
        end
    end

    % Interpolate Directivity
    nInterp = 360; % number of interpolation points
    interpFactor = max(round(nInterp/length(sigma)),1);
    interpPoints = (1:1/interpFactor:length(sigma))';
    sigma_int = interp1(sigma,interpPoints,'linear');
    directivity_int = interp1(directivity_norm,interpPoints,'pchip');

    % PLOT VERTICAL DIRECTIVITY
    sigma_plo = [sigma_int; sigma_int(1)];
    directivity_plo = [directivity_int; directivity_int(1,:)];
    h = polarPlot(sigma_plo,directivity_plo,[0 1],legendText);
    title({'Horizontal Directivity of Air Gun Array';...
        sprintf('\\rm\\theta = %0.0f',below360(targetElevation,'deg'));'';''})
end

