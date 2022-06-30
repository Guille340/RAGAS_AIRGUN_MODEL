% Radiation = bubbleRadiation(Bubble,receiverPositions,varargin)
%
% DESCRIPTION
% Computes the acoustic signature of the air gun array at the specified 
% RECEIVERPOSITIONS (pressure and particle velocity), and related metrics 
% (broadband pressure RMS/SEL/Peak/Peak-to-Peak, band pressure RMS/SEL, 
% primary and bubble amplitudes and times). The results and configuration 
% information are returned as a RADIATION structure.
%
% BUBBLERADIATION.m calls the functions RADIATIONFIRSTORDER.m or 
% RADIATIONSECONDORDER.m to compute the time-dependent pressure and particle
% velocity using the first-order or second-order approximations from
% Gilmore (1952). The metrics are calculated with SIGNATUREMETRICS.m.
%
% INPUT ARGUMENTS
% - Bubble: structure containing the input configuration information and 
%   predicted time-depended bubble motion parameters. The structure is 
%   generated with functions MULTIPHYSICSDYNAMICS.m or POLYTROPICDYNAMICS.m.
% - receiverPositions: 3-element vector or matrix of dimensions [3,nReceivers]
%   containing the [x y z] position of the receiver(s).
%
% INPUT PROPERTIES
% - 'GunIdentifiers': identifiers of the air guns in the array to be 
%    considered for the processing of the acoustic signature. By default,
%    BUBBLERADIATION.m works with the air guns selected in BUBBLE.CONFIG.
%    The 'GunIndentifiers' property is useful, for example, when only one air 
%    gun is to be processed (e.g., when computing bubble interactions as a
%    result of the acoustic emission of individual air guns).
% - 'Ghost': TRUE for the ghost reflection to be considered for the acoustic
%    radiation.
% - 'DisplayProgress': TRUE for displaying the progress of the simulations in
%    the command window.
%
% OUTPUT ARGUMENTS
% - Radiation: structure containing the acoustic signature of the array at
%   the specified receiver location, related metrics, and configuration 
%   information. For details about the content of this structure see funtion
%   INITIALISERADIATIONSTRUCT.m
%
% FUNCTION CALL
% 1. Radiation = bubbleRadiation(Bubble,receiverPositions)
% 2. Radiation = bubbleRadiation(...,PROPERTY,VALUE)
%    Properties: 'GunIdentifiers', 'Ghost', 'DisplayProgress'
% 
% FUNCTION DEPENDENCIES
% - bubbleRadiation_ec
% - arrayCentre
% - radiationFirstOrder
% - radiationSecondOrder
% - signatureMetrics
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% REFERENCES
% Gilmore, F. R. (1952) "The growth or collapse of a spherical bubble in a 
% viscous compressible liquid", California Institute of Technology, report 26-4
%
% Ziolkowski, A. (1970). "A method for calculating the output pressure 
% waveform from an air gun”, Geophys. J. R. Astr. Soc. 21(2), 137-161
%
% See also RADIATIONFIRSTORDER, RADIATIONSECONDORDER, SIGNATUREMETRICS

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function Radiation = bubbleRadiation(Bubble,receiverPositions,varargin)

narginchk(2,8)

% ERROR CONTROL
[gunIds,ghost,displayProgress] = bubbleRadiation_ec(Bubble,...
    receiverPositions,varargin);

% Reshape 'receiverPositions'
[M,N] = size(receiverPositions);
if M ~= N && (find([M,N] == 3) == 2)
    receiverPositions = receiverPositions';
end

% CONFIGURATION PARAMETERS
tmax = Bubble.Config.duration;
fs = Bubble.Config.sampleRate;
order = 2; % order of radiation formula

% GENERAL VARIABLES
[x_arr,y_arr,z_arr] = arrayCentre(Bubble); % central position of air gun array
nPoints = tmax*fs + 1; % number of points in the airgun signature
nGuns = length(gunIds);
nReceivers = size(receiverPositions,2);

% DISPLAY PROGRESS 
if displayProgress
    fprintf('\nPROCESSING ACOUSTIC RADIATION\n')
end

for m = 1:nReceivers
    
    % Display Progress
    if displayProgress
        fprintf('# Pressure and velocity at receiver (%d/%d)\n',m,nReceivers)
    end
    
    % Receiver Location
    x_rec = receiverPositions(1,m);
    y_rec = receiverPositions(2,m);
    z_rec = receiverPositions(3,m);
    
    % Source to Receiver Distance and Angles
    dist = sqrt((x_rec - x_arr)^2 + (y_rec - y_arr)^2 + (z_rec - z_arr)^2); % distance [m]
    azim = atan2((x_rec - x_arr),(y_rec - y_arr)) * 180/pi; % azimuth [deg]
    elev = asin((z_rec - z_arr)/dist) * 180/pi; % elevation [deg]

    % Individual Notional Signatures
    for q = 1:nGuns
        % Source Location
        iGun = find([Bubble.Array.gunId] == gunIds(q));
        x_sou = Bubble.Array(iGun).gunX;
        y_sou = Bubble.Array(iGun).gunY;
        z_sou = Bubble.Array(iGun).gunZ;

        % Source-Receiver Distances (direct and ghost)
        Rd = sqrt((x_rec - x_sou)^2 + (y_rec - y_sou)^2 + (z_rec - z_sou)^2);
        Rr = sqrt((x_rec - x_sou)^2 + (y_rec - y_sou)^2 + (z_rec + z_sou)^2);

        % Pressure and Particle Velocity
        if order == 1
            [pd(:,q),ud(:,q),td(1,q)] = radiationFirstOrder(Bubble.Dynamics(iGun),Rd);
            [pr(:,q),ur(:,q),tr(1,q)] = radiationFirstOrder(Bubble.Dynamics(iGun),Rr);   
        else % order == 2
            [pd(:,q),ud(:,q),td(1,q)] = radiationSecondOrder(Bubble.Dynamics(iGun),Rd,fs);
            [pr(:,q),ur(:,q),tr(1,q)] = radiationSecondOrder(Bubble.Dynamics(iGun),Rr,fs);   
        end
    end

    % Combined Notional Signature
    tmin = min(td);
    td = td - tmin;
    tr = tr - tmin;
    p = zeros(nPoints,1);
    u = zeros(nPoints,1);
    for q = 1:nGuns
        if ghost
            nzd = round(td(q)*fs); % number of shifting zeroes (direct)
            nzr = round(tr(q)*fs); % number of shifting zeroes (reflected)
            pdz = [zeros(nzd,1); pd(1:end - nzd,q)]; 
            prz = [zeros(nzr,1); pr(1:end - nzr,q)];
            udz = [zeros(nzd,1); ud(1:end - nzd,q)]; 
            urz = [zeros(nzr,1); ur(1:end - nzr,q)]; 
            p = p + pdz - prz;
            u = u + udz - urz;
        else 
            nzd = round(td(q)*fs); % number of shifting zeroes (direct)
            pdz = [zeros(nzd,1) ; pd(1:end - nzd,q)]; 
            udz = [zeros(nzd,1) ; ud(1:end - nzd,q)]; 
            p = p + pdz;
            u = u + udz;
        end
    end
    
    % Normalised Pressure and Velocity Signatures (ref. array centre)
    p = p.*dist; % pressure normalised to 1 m
    u = u.*dist; % particle velocity normalised to 1 m
    
    % Acoustic Parameters
    freqLimits = [10 fs/2];
    Metrics = signatureMetrics(p,fs,freqLimits);
    
    % Populate 'Radiation' Structure
    Radiation.Config = Bubble.Config;

    Radiation.Array = Bubble.Array;
    
    Radiation.Signature(m).timeDelay = tmin; % alternatively, delay from array centre is dist/cw
    Radiation.Signature(m).pressureAtOneMetre = p;
    Radiation.Signature(m).velocityAtOneMetre = u;
    Radiation.Signature(m).receiverPosition = [x_rec, y_rec, z_rec];
    Radiation.Signature(m).sourcePosition = [x_arr, y_arr, z_arr];
    Radiation.Signature(m).souToRecDistance = round(dist*1e3)*1e-3; % distance (round to 1 mm)
    Radiation.Signature(m).souToRecAzimuth = round(azim*1e3)*1e-3; % azimuth (round to 1e-3 deg)
    Radiation.Signature(m).souToRecElevation = round(elev*1e3)*1e-3; % elevation (round to 1e-3 deg)
    
    Radiation.Metrics(m) = Metrics;
end
        