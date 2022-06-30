% Plot = readPlot(Root)
%
% DESCRIPTION
% Reads and verifies the configuration information from the "plot.json" file
% stored in <ROOT.DATAFOLDER> and returns that information as a PLOT 
% structure. 
%
% "plot.json" is a JSON structure. The values need to be introduced manually
% in accordance with the JSON format. The file contains 7 fields necessary
% for the computation and representation of the acoustic signatures,
% array directivity and spectrum. Below is an example of "plot.json".
%
%  {
%  "souToRecDistance": 1000,
%  "souToRecAzimuth": 0,
%  "souToRecElevation": 20,
%  "frequency": [50,100,200],
%  "radiationFormula": 2,
%  "interaction": false,
%  "radiationGhost": true
% }
% 
% INPUT ARGUMENTS
% - Root: structure containing root directory information. Currently, ROOT
%   only contains folder where the configuration files and output data is
%   stored (ROOT.DATAFOLDER).
%
% OUTPUT ARGUMENTS
% - Plot: single-element structure containing the configuration information
%   necessary to compute and represent the acoustic signature, directivity
%   and spectrum of the array. For details about the contents of this structure
%   see function INITIALISEPLOTSTRUCT.m.
%
% DEPENDENCIES
% - None
%
% FUNCTION CALL
% 1. Plot = readPlot(Root)
% 
% FUNCTION DEPENDENCIES
% - None
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% See also READCONFIG

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function Plot = readPlot(Root)

narginchk(1,1) % check number of input arguments

% INITIALISE CONFIGURATION STRUCTURE
Plot = initialisePlotStruct();
Plot.inputStatus = true;

% CHECK IF ROOT FOLDER EXISTS
if exist(Root.dataFolder,'dir') == 7
    plotPath = fullfile(Root.dataFolder,'plot.json');
    
    % Read Plot File
    if exist(plotPath,'file') == 2
        PlotFile = jsondecode(fileread(plotPath)); % read settings
    else
        Plot.inputStatus = false;
        warning(['The plot file PLOT.JSON could not be found '...
            'in folder %s'],Root.dataFolder)
    end
else
    PlotFile = [];
    Plot.inputStatus = false;
    warning('The root folder in ROOT.JSON does not exist')
end

% VERIFY INDIVIDUAL FIELDS IN 'plot.json'
if isstruct(PlotFile)
    fieldNames = fieldnames(PlotFile);
    nFields = numel(fieldNames); % number of fields in temporal structure
    for m = 1:nFields
        fieldName = fieldNames{m}; % current field name
        fieldValue = [PlotFile.(fieldName)]; % current field value

        switch fieldName
            case 'souToRecDistance'
                if ~isnumeric(fieldValue) || ~isscalar(fieldValue) ...
                        || fieldValue < 0
                    fieldValue = 1000;
                    warning(['SOUTORECDISTANCE must be a positive scalar. '...
                        'SOUTORECDISTANCE = %0.1f will be used'],fieldValue)
                end
            case 'souToRecAzimuth'
                if ~isnumeric(fieldValue) || ~isscalar(fieldValue)
                    fieldValue = 0;
                    warning(['SOUTORECAZIMUTH must be a positive scalar. '...
                        'SOUTORECAZIMUTH = %0.1f will be used'],fieldValue)
                elseif fieldValue < -180 || fieldValue > 360 %#ok<*BDSCI>
                    fieldValue = below360(fieldValue,'deg');
                    warning(['SOUTORECAZIMUTH must be within the ranges '...
                        '-180 to 180 deg or 0 to 360 deg. SOUTORECAZIMUTH '...
                        '= %0.0f will be used'],fieldValue)
                end
            case 'souToRecElevation'
                if ~isnumeric(fieldValue) || ~isscalar(fieldValue)
                    fieldValue = 45;
                    warning(['SOUTORECELEVATION must be a positive scalar. '...
                        'SOUTORECELEVATION = %0.1f will be used'],fieldValue)
                elseif fieldValue < 0 || fieldValue > 90 %#ok<*BDSCI>
                    fieldValue = below360(fieldValue,'deg');
                    warning(['SOUTORECELEVATION must be within the range '...
                        '0 to 90 deg. SOUTORECELEVATION = %0.0f will be '...
                        'used'],fieldValue)
                end
            case 'frequency'
                if ~isnumeric(fieldValue) || ~isvector(fieldValue) ...
                        || any(fieldValue < 0)
                    fieldValue = [];
                    warning(['FREQUENCY must be a numeric positive vector. '...
                        'FREQUENCY = [] will be used (broadband)'])
                end  
            case 'radiationFormula'
                if ~isnumeric(fieldValue) || ~isscalar(fieldValue) ...
                        || any(fieldValue ~= [1 2])
                    fieldValue = 2;
                    warning(['RADIATIONFORMULA must be 1 or 2. '...
                        'RADIATIONFORMULA = ''%d'' will be used'],fieldValue)
                end
            case 'interaction'
                if isempty(fieldValue) || (~isnumeric(fieldValue) ...
                    && ~islogical(fieldValue)) || ~any(fieldValue == [0 1])
                    fieldValue = true;
                    warning(['INTERACTION must be [0 1] or logical. '...
                        'INTERACTION = TRUE will be used'])
                end 
            case 'radiationGhost'
                if isempty(fieldValue) || ~any(fieldValue == [0 1])
                    fieldValue = true;
                    warning(['RADIATIONGHOST must be [0 1] or logical. '...
                        'RADIATIONGHOST = TRUE will be used'])
                end 
        end 
        Plot(1).(fieldName) = fieldValue;
    end
else
    warning('PLOTFILE is not a valid plotting structure')
end
