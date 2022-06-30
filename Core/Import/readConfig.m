% Config = readConfig(Root)
%
% DESCRIPTION
% Reads and verifies the configuration information from the "config.json" file
% stored in <ROOT.DATAFOLDER> and returns that information as a CONFIG 
% structure. The function also calls READARRAY.m to extract the layout of
% the air gun array and incorporates that information into the CONFIG structure.
%
% "config.json" is a JSON structure. The values need to be introduced manually
% in accordance with the JSON format. The file contains 10 fields necessary
% for the computation of the bubble dynamics of the array. Below is an example
% of "config.json".
%
%  {
%  "gunIds": [1,2,3,4],
%  "duration": 0.5,
%  "sampleRate": 10000,
%  "waterDensity": [],
%  "waterSoundSpeed": [],
%  "waterTemperature": [],
%  "waterSalinity": [],
%  "bubbleMethod": "multiphysics",
%  "bubbleFormula": "gilmore",
%  "interactionGhost": false
% 
% INPUT ARGUMENTS
% - Root: structure containing root directory information. Currently, ROOT
%   only contains folder where the configuration files and output data is
%   stored (ROOT.DATAFOLDER).
%
% OUTPUT ARGUMENTS
% - Config: single-element structure containing the configuration information
%   necessary to process the bubble dynamics of the array. For details about 
%   the contents of this structure see function INITIALISECONFIGSTRUCT.m.
%
% DEPENDENCIES
% - None
%
% FUNCTION CALL
% 1. Config = readConfig(Root)
% 
% FUNCTION DEPENDENCIES
% - readArray
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% See also READARRAY, READPLOT

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function Config = readConfig(Root)

narginchk(1,1) % check number of input arguments

% INITIALISE CONFIGURATION STRUCTURE
Config = initialiseConfigStruct();
Config.inputStatus = true;

% CHECK IF ROOT FOLDER EXISTS
if exist(Root.dataFolder,'dir') == 7
    configPath = fullfile(Root.dataFolder,'config.json');
    
    % Read Configuration File
    if exist(configPath,'file') == 2
        ConfigFile = jsondecode(fileread(configPath)); % read settings
    else
        Config.inputStatus = false;
        warning(['The configuration file CONFIG.JSON could not be found '...
            'in folder %s'],Root.dataFolder)
    end
    
    % Read Air Gun Array Layout Data
    Config.Array = readArray(Root);
else
    ConfigFile = [];
    Config.inputStatus = false;
    warning('The root folder in ROOT.JSON does not exist')
end

% GUN IDS
gunIds = [];
if ~isempty(Config.Array)
    gunIds = [Config.Array.gunId];
end

% VERIFY INDIVIDUAL FIELDS IN 'config.json'
if isstruct(ConfigFile)
    fieldNames = fieldnames(ConfigFile);
    nFields = numel(fieldNames); % number of fields in temporal structure
    for m = 1:nFields
        fieldName = fieldNames{m}; % current field name
        fieldValue = [ConfigFile.(fieldName)]; % current field value

        switch fieldName
            case 'gunIds'
                if ~isempty(fieldValue)
                    if ~isnumeric(fieldValue) || ~isvector(fieldValue) ...
                            || any(rem(fieldValue,1))
                        fieldValue = gunIds;
                        warning(['GUNIDS must be a vector of integers. '...
                            'All air guns in the array will be used'])
                    elseif ~isempty(gunIds) && any(~ismember(fieldValue,gunIds))
                        fieldValue = fieldValue(ismember(fieldValue,gunIds));
                        fieldValue = fieldValue(:)';
                        warning(['One or more values in GUNIDS is '...
                            'not included in the array. Only those gun ids '...
                            'in the array will be considered'])
                    else
                        fieldValue = fieldValue(:)';
                    end
                else
                    fieldValue = gunIds;
                end  
            case 'duration'
                if ~isnumeric(fieldValue) || ~isscalar(fieldValue) ...
                        || fieldValue < 0 %#ok<*BDSCI>
                    fieldValue = 1;
                    warning(['DURATION must be a positive scalar. '...
                        'DURATION = %0.1f will be used'],fieldValue)
                end 
            case 'sampleRate'
                if ~isnumeric(fieldValue) || ~isscalar(fieldValue) ...
                        || fieldValue < 0
                    fieldValue = 1e4;
                    warning(['SAMPLERATE must be a positive scalar. '...
                        'SAMPLERATE = %0.1f will be used'],fieldValue)
                elseif fieldValue < 200
                    warning(['SAMPLERATE must be larger than 200 Hz '...
                        'for minimally accurate results'])
                end  
            case 'waterDensity'
                if ~isnumeric(fieldValue) || ~isscalar(fieldValue) ...
                        || fieldValue < 0
                    fieldValue = 1024;
                    warning(['WATERDENSITY must be a positive scalar. '...
                        'WATERDENSITY = %0.1f will be used'],fieldValue)
                end 
            case 'waterSoundSpeed'
                if ~isnumeric(fieldValue) || ~isscalar(fieldValue) ...
                        || fieldValue < 0
                    fieldValue = 1489;
                    warning(['WATERSOUNDSPEED must be a positive scalar. '...
                        'WATERSOUNDSPEED = %0.1f will be used'],fieldValue)
                end  
            case 'waterTemperature'
                if ~isnumeric(fieldValue) || ~isscalar(fieldValue) ...
                        || fieldValue < 0
                    fieldValue = 284.6;
                    warning(['WATERTEMPERATURE must be a positive scalar. '...
                        'WATERTEMPERATURE = %0.1f will be used'],fieldValue)
                end
            case 'waterSalinity'
                if ~isnumeric(fieldValue) || ~isscalar(fieldValue) ...
                        || fieldValue < 0
                    fieldValue = 35;
                    warning(['WATERSALINITY must be a positive scalar. '...
                        'WATERSALINITY = %0.1f will be used'],fieldValue)
                end 
            case 'bubbleMethod'
                if ~ischar(fieldValue) || ~ismember(lower(fieldValue),...
                        {'multiphysics','polytropic'})
                    fieldValue = 'multiphysics';
                    warning(['The option introduced for BUBBLEMETHOD '...
                        'is not valid. BUBBLEMETHOD = ''%s'' '...
                        'will be used'],fieldValue);
                end
            case 'bubbleFormula'
                if ~ischar(fieldValue) || ~ismember(lower(fieldValue),...
                        {'gilmore','keller-kollodner','rayleigh-plesset',...
                        'herring'})
                    fieldValue = 'gilmore';
                    warning(['The option introduced for BUBBLEFORMULA '...
                        'is not valid. BUBBLEFORMULA = ''%s'' '...
                        'will be used'],fieldValue)
                end
            case 'interactionGhost'
                if isempty(fieldValue) || ~any(fieldValue == [0 1])
                    fieldValue = true;
                    warning(['INTERACTIONGHOST must be [0 1] or logical. '...
                        'INTERACTIONGHOST = TRUE will be used'])
                end
        end
        Config(1).(fieldName) = fieldValue;
    end
else
    warning('CONFIGFILE is not a valid configuration structure')
end
