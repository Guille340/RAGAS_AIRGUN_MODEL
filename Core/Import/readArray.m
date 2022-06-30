% Array = readArray(Root)
%
% DESCRIPTION
% Reads and verifies the air gun array layout data from the "array.csv"
% table stored in <ROOT.DATAFOLDER> and returns that information as an ARRAY 
% structure.
%
% "array.csv" is a table with comma-separated values consisting in a first-line
% header followed by one data line per air gun in the array. Below is an 
% example of an array table. The order of the fields is irrelevant. The field
% names in the header are case-insensitive (upper/lower).
% 
% Id,Pressure_psi,GeneratorVolume_cuin,InjectorVolume_cuin,...
%    InjectorTime_ms,Radius_m,X_m,Y_m,Z_m,
% 1,2000.0,45.0,0.0,0.0,0.2,0.000,0.500,6.000, 
% 2,2000.0,45.0,0.0,0.0,0.2,0.000,-0.500,6.000,
% 3,2000.0,70.0,0.0,0.0,0.2,1.700,0.500,6.000,
% 4,2000.0,70.0,0.0,0.0,0.2,1.700,-0.500,6.000,
%
% INPUT ARGUMENTS
% - Root: structure containing root directory information. Currently, ROOT
%   only contains folder where the configuration files and output data is
%   stored (ROOT.DATAFOLDER).
%
% OUTPUT ARGUMENTS
% - Array: multi-element structure containing the properties of air guns in
%   the seismic array. Each element corresponds to one air gun in the array.
%   For details about the contents of this structure see function
%   INITIALISEARRAYSTRUCT.m.
%
% DEPENDENCIES
% - None
%
% FUNCTION CALL
% 1. Array = readArray(Root)
% 
% FUNCTION DEPENDENCIES
% - None
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% See also READCONFIG, READPLOT

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function Array = readArray(Root)

% Check if Root Folder Exists
if exist(Root.dataFolder,'dir') == 7
    arrayPath = fullfile(Root.dataFolder,'array.csv');
       
    % Read Air Gun Array Layout Data
    if exist(arrayPath,'file') ~= 2
        warning(['The array layout file ARRAY.CSV could not be found '...
            'in folder %s'],Root.dataFolder)
    end
else
    warning('The root folder in ROOT.JSON does not exist')
end

% Read Array Layout File
Array = [];
fid = fopen(arrayPath,'r');
if fid ~= -1
    % Read Data from 'array.csv'
    header = textscan(fid,'%s%s%s%s%s%s%s%s%s',1,'Delimiter',',');
    data = textscan(fid,'%s%s%s%s%s%s%s%s%s','Delimiter',',','HeaderLines',1);
    nGuns = length(data{1});
    nFields = length(data);
    
    % Check if File is an Array Layout File
    fieldNames_valid = {'id','pressure_psi','generatorvolume_cuin',...
        'injectorvolume_cuin','injectortime_ms','radius_m','x_m','y_m','z_m'};
    fieldNames = repmat({''},1,nFields);
    for m = 1:nFields
        fieldNames{m} = lower(header{m}{1});
    end
    if any(~ismember(fieldNames,fieldNames_valid)) ...
            || any(~ismember(fieldNames_valid,fieldNames))
        error('''array.csv'' is not a valid array layout structure')        
    end
    
    % Populate Array Structure
    Array = initialiseArrayStruct();   
    for m = 1:nFields
        switch fieldNames{m}
            case 'id'
                % Check Validity
                fieldData = data{m};
                isValid = false(1,nGuns);
                for n = 1:nGuns
                    id = str2num(fieldData{n}); %#ok<*ST2NM>
                    isValid(n) = ~isempty(id) && isnumeric(id) ...
                        && isscalar(id) && ~rem(id,1);                        
                end
                if any(~isValid)
                    warning(['The value of the ''Id'' field in '...
                        '''array.csv'' for one or more air guns is not '...
                        'an integer scalar'])
                end  
                
                % Populate 'id' field of Array Structure
                for n = 1:nGuns
                    if isValid(n)
                        Array(n).gunId = str2num(fieldData{n});
                    end
                end   
                
            case 'pressure_psi'
                % Check Validity
                fieldData = data{m};
                isValid = false(1,nGuns);
                for n = 1:nGuns
                    pressure = str2num(fieldData{n});
                    isValid(n) = ~isempty(pressure) && isnumeric(pressure) ...
                        && isscalar(pressure) && pressure > 0;
                end
                if any(~isValid)
                    warning(['The value of the ''Pressure'' field in '...
                        '''array.csv'' for one or more air guns is not '...
                        'a positive scalar'])
                end  
                
                % Populate 'pressure' field of Array Structure
                for n = 1:nGuns
                    if isValid(n)
                        Array(n).gunPressure = str2num(fieldData{n});
                    end
                end 
                
            case 'generatorvolume_cuin'
                % Check Validity
                fieldData = data{m};
                isValid = false(1,nGuns);
                for n = 1:nGuns
                    volumeGenerator = str2num(fieldData{n});
                    isValid(n) = ~isempty(volumeGenerator) ...
                        && isnumeric(volumeGenerator) ...
                        && isscalar(volumeGenerator) && volumeGenerator > 0;
                end
                if any(~isValid)
                    warning(['The value of the ''VolumeGenerator'' field '...
                        'in ''array.csv'' for one or more air guns is not '...
                        'a positive scalar'])
                end  
                
                % Populate 'volumeGenerator' field of Array Structure
                for n = 1:nGuns
                    if isValid(n)
                        Array(n).gunGeneratorVolume = str2num(fieldData{n});
                    end
                end

            case 'injectorvolume_cuin'
                % Check Validity
                fieldData = data{m};
                isValid = false(1,nGuns);
                for n = 1:nGuns
                    volumeInjector = str2num(fieldData{n});
                    isValid(n) = ~isempty(volumeInjector) ...
                        && isnumeric(volumeInjector) ...
                        && isscalar(volumeInjector) && volumeInjector >= 0;
                end
                if any(~isValid)
                    warning(['The value of the ''VolumeInjector'' field '...
                        'in ''array.csv'' for one or more air guns is not '...
                        'a positive scalar or zero'])
                end  
                
                % Populate 'volumeInjector' field of Array Structure
                for n = 1:nGuns
                    if isValid(n)
                        Array(n).gunInjectorVolume = str2num(fieldData{n});
                    end
                end
                
            case 'injectortime_ms'
                % Check Validity
                fieldData = data{m};
                isValid = false(1,nGuns);
                for n = 1:nGuns
                    timeInjector = str2num(fieldData{n}); 
                    isValid(n) = ~isempty(timeInjector) ...
                        && isnumeric(timeInjector) ...
                        && isscalar(timeInjector) && timeInjector >= 0;
                end
                if any(~isValid)
                    warning(['The value of the ''TimeInjector'' field '...
                        'in ''array.csv'' for one or more air guns is not '...
                        'a positive scalar or zero'])
                end  
                
                % Populate 'timeInjector' field of Array Structure
                for n = 1:nGuns
                    if isValid(n)
                        Array(n).gunInjectorTime = str2num(fieldData{n});
                    end
                end

            case 'radius_m'
                % Check Validity
                fieldData = data{m};
                isValid = false(1,nGuns);
                for n = 1:nGuns
                    radius = str2num(fieldData{n});
                    isValid(n) = ~isempty(radius) && isnumeric(radius) ...
                        && isscalar(radius) && radius >= 0;
                end
                if any(~isValid)
                    warning(['The value of the ''Radius'' field '...
                        'in ''array.csv'' for one or more air guns is not '...
                        'a positive scalar or zero'])
                end  
                
                % Populate 'radius' field of Array Structure
                for n = 1:nGuns
                    if isValid(n)
                        Array(n).gunRadius = str2num(fieldData{n});
                    end
                end
                
            case 'x_m'
                % Check Validity
                fieldData = data{m};
                isValid = false(1,nGuns);
                for n = 1:nGuns
                    x = str2num(fieldData{n});
                    isValid(n) = ~isempty(x) && isnumeric(x) && isscalar(x);
                end
                if any(~isValid)
                    warning(['The value of the ''X'' field in '...
                        '''array.csv'' for one or more air guns is '...
                        'not a scalar'])
                end  
                
                % Populate 'x' field of Array Structure
                for n = 1:nGuns
                    if isValid(n)
                        Array(n).gunX = str2num(fieldData{n});
                    end
                end
                
            case 'y_m'
                % Check Validity
                fieldData = data{m};
                isValid = false(1,nGuns);
                for n = 1:nGuns
                    y = str2num(fieldData{n}); %#ok<*ST2NM>
                    isValid(n) = ~isempty(y) && isnumeric(y) ...
                        && isscalar(y);
                end
                if any(~isValid)
                    warning(['The value of the ''Y'' field in '...
                        '''array.csv'' for one or more air guns is '...
                        'not a scalar'])
                end  
                
                % Populate 'y' field of Array Structure
                for n = 1:nGuns
                    if isValid(n)
                        Array(n).gunY = str2num(fieldData{n});
                    end
                end

            case 'z_m'
                % Check Validity
                fieldData = data{m};
                isValid = false(1,nGuns);
                for n = 1:nGuns
                    z = str2num(fieldData{n}); %#ok<*ST2NM>
                    isValid(n) = ~isempty(z) && isnumeric(z) && isscalar(z);
                end
                if any(~isValid)
                    warning(['The value of the ''Z'' field in '...
                        '''array.csv'' for one or more air guns is '...
                        'not a scalar'])
                end  
                
                % Populate 'z' field of Array Structure
                for n = 1:nGuns
                    if isValid(n)
                        Array(n).gunZ = str2num(fieldData{n});
                    end
                end
        end 
    end 
end
fclose(fid);

% Update Input Status
for n = 1:nGuns
    Array(n).inputStatus = ~isempty(Array(n).gunId) ...
        && ~isempty(Array(n).gunPressure) ...
        && ~isempty(Array(n).gunGeneratorVolume) ...
        && ~isempty(Array(n).gunInjectorVolume) ...
        && ~isempty(Array(n).gunInjectorTime) ...
        && ~isempty(Array(n).gunRadius) ...
        && ~isempty(Array(n).gunX) ...
        && ~isempty(Array(n).gunY) ...
        && ~isempty(Array(n).gunZ); %#ok<AGROW>
end
