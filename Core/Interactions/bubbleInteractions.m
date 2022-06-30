% Interactions = bubbleInteractions(Config,Bubble,varargin)
%
% DESCRIPTION
% Computes the effective sound pressure produced by the surrounding air guns
% at every air gun location and returns that information as an INTERACTIONS
% structure. The effective pressure can be later added to the hydrostatic 
% pressure to account for the effect that nearby bubbles will have on the 
% motion of each individual bubble, as described by Ziolkowsky et al. (1982).
% BUBBLEINTERACTIONS.m calls function BUBBLERADIATION.m to calculate the
% acoustic signature produced by each individual air gun acting independently
% on each bubble.
%
% The functions MULTIPHYSICSDYNAMICS.m and POLYTROPICDYNAMICS.m accept an
% optional INTERACTIONS structure to simulate the bubble motion including
% interactions.
%
% INPUT ARGUMENTS
% - Config: configuration structure generated from the JSON configuration
%   file stored in <ROOT.DATAFOLDER> directory (config.json). For details
%   about the content of this structure see INITIALISECONFIGSTRUCT.m.
% - Bubble: structure containing the input configuration information and 
%   predicted time-depended bubble motion parameters. The structure is 
%   generated with functions MULTIPHYSICSDYNAMICS.m or POLYTROPICDYNAMICS.m.
%
% INPUT PROPERTIES
% - 'DisplayProgress': TRUE for displaying the progress of the simulations in
%    the command window.
%
% OUTPUT ARGUMENTS
% - Interactions: structure containing the effective acoustic pressure 
%   produced on each bubble by the surrounding nearby bubbles. The structure
%   also includes the identifiers of the various interacting guns. For details
%   about the content of this structure see INITIALISEINTERACTIONSSTRUCT.m
%
% FUNCTION CALL
% 1. Interactions = bubbleInteractions(Config,Bubble)
% 2. Interactions = bubbleInteractions(...,PROPERTY,VALUE)
%    Properties: 'DisplayProgress'
% 
% FUNCTION DEPENDENCIES
% - bubbleInteractions_ec
% - initialiseInteractionsStruct
% - bubbleRadiation
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% REFERENCES
% Ziolkowski,, A., Parkes, G., Hatton, L. and Haugland, T. (1982). "The 
% signature of an air gun array: computation from near-field measurements 
% including interactions", Geophysics 47(10), 1413-1421
%
% See also INITIALISEINTERACTIONSSTRUCT, BUBBLERADIATION,
% MULTIPHYSICSDYNAMICS, POLYTROPICDYNAMICS

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function Interactions = bubbleInteractions(Config,Bubble,varargin)

narginchk(2,4)

% ERROR CONTROL
displayProgress = bubbleInteractions_ec(Config,Bubble,varargin);

% GENERAL PARAMETERS
gunIds = Config.gunIds;
interactionGhost = Config.interactionGhost;
tmax = Config.duration; % signature duration [s]
fs = Config.sampleRate; % sampling frequency for returned parameters [Hz]

% DISPLAY PROGRESS 
if displayProgress
    fprintf('\nPROCESSING BUBBLE INTERACTIONS\n')
end

% CHECK COALESCENSE
if displayProgress
    fprintf('# Checking coalescense\n')
end
if isCoalescent(Bubble)
    warning(['One or more bubble pairs experience coalescence. The '...
        'accuracy of the results may be limited'])
end

% PROCESS EFFECTIVE SOUND PRESSURE AT EACH AIR GUN LOCATION
nPoints = round(tmax * fs + 1);
nGuns = length(gunIds);
Interactions = initialiseInteractionsStruct(Config);
for q = 1:nGuns
    % Display Progress
    if displayProgress
        fprintf('# Processing radiation \n')
    end

    % Load Parameters
    gunId_rec = gunIds(q); % id of airgun affected by surrounding airguns
    gunIds_sou = gunIds(gunIds ~= gunId_rec); % id of airguns interacting with current airgun

    % Calculate Position of Receiving Airgun
    iGun_rec = q;
    x_rec = Bubble.Array(iGun_rec).gunX; % horizontal shift from the centre of the array (receiving airgun) [m]
    y_rec = Bubble.Array(iGun_rec).gunY; % vertical shift from the centre of the array (receiving airgun) [m]
    z_rec = Bubble.Array(iGun_rec).gunZ; % depth (receiving airgun) [m]

    % Calculate Sound Pressure at Receiving Airgun from Surrounding Ones
    tmin_sou = 0;
    ps_sou = zeros(nPoints,1);
    if ~isempty(gunIds_sou)
        Radiation = bubbleRadiation(Bubble,[x_rec y_rec z_rec],...
            'GunIdentifiers',gunIds_sou,'Ghost',interactionGhost,...
            'DisplayProgress',false); %#ok<*FNDSB>
        ps_sou = Radiation.Signature.pressureAtOneMetre ...
            /Radiation.Signature.souToRecDistance;
        tmin_sou = Radiation.Signature.timeDelay;
    end
    
    % Calculate Sound Pressure at Receiving Airgun from Itself (Ghost)
    tmin_rec = 0;
    ps_rec = zeros(nPoints,1);
    if interactionGhost
        Radiation = bubbleRadiation(Bubble,[x_rec y_rec -z_rec],...
            'GunIdentifiers',gunId_rec,'Ghost',false,...
            'DisplayProgress',false); %#ok<*FNDSB>
        ps_rec = Radiation.Signature.pressureAtOneMetre ...
            /Radiation.Signature.souToRecDistance;
        tmin_rec = Radiation.Signature.timeDelay;
    end
    
    % Calculate Sound Pressure at Receiver Airgun from All Contributions
    if ~tmin_sou && ~tmin_rec
        tmin = min([tmin_sou tmin_rec]);
        tmin_sou = tmin_sou - tmin;
        tmin_rec = tmin_rec - tmin;
    else % if any delay is zero, set both as zero
        tmin_sou = 0;
        tmin_rec = 0;
    end
    nz_sou = round(tmin_sou*fs); % number of shifting zeroes (direct)
    nz_rec = round(tmin_rec*fs); % number of shifting zeroes (reflected)
    ps_sou = [zeros(nz_sou,1); ps_sou(1:end - nz_sou)]; 
    ps_rec = [zeros(nz_rec,1); ps_rec(1:end - nz_rec)];
    ps = ps_sou - ps_rec;

    % Build Structure of Effective Sound Pressure
    Interactions(q).gunId = gunId_rec;
    Interactions(q).gunIds_sou = gunIds_sou;
    Interactions(q).gunIds_all = gunIds;
    Interactions(q).dynamicPressure = ps;
end

