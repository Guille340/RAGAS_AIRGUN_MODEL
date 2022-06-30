function [Interactions,Calibration,excludePhysics,displayProgress] ...
    = polytropicDynamics_ec(Config,varargin_cell)

% Verify number of Input Arguments
nFixArg = 1;
varargin = varargin_cell;
nVarargin = length(varargin);
nargin = nFixArg + nVarargin;
if rem(nargin - nFixArg,2)
    error('Variable input arguments must come in pairs (PROPERTY,VALUE)')
end

% Error Control ('Config')
if ~isConfigStruct(Config)
    error('Input CONFIG is not a valid configuration structure')
end

% Parameters from 'Config' Structure
tmax = Config.duration; % duration of parameter waveforms [s]
fs = Config.sampleRate; % sampling rate of bubble parameters [Hz]
fsp = round(1e4/fs)*fs; % sampling rate for parameter processing [Hz]
gunIds = Config.gunIds; % identifiers of selected air guns
nPoints = round(tmax*fsp + 1); % number of points in the airgun signature

% Extract and Verify Input Properties
validProperties = {'interactions','calibration','excludephysics',...
    'displayprogress'};
properties = lower(varargin(1:2:end));
if any(~ismember(properties,validProperties))
    error('One or more input properties are not recognised')
end

% Default Input Values
Interactions = initialiseInteractionsStruct(Config);
Calibration = [];
excludePhysics = '';
displayProgress = false;

% Extract and Verify Input Values
values = varargin(2:2:end);
nPairs = (nargin - nFixArg)/2; % number of (PROPERTY,VALUE) pairs
for m = 1:nPairs
    property = properties{m};
    switch property % populate with more properties if needed
        case 'interactions' % effective hydrostatic pressure
            Interactions = values{m};
            if isInteractionsStruct(Interactions)
                if any(~ismember(gunIds,[Interactions.gunId])) ...
                        || any(~ismember([Interactions.gunId],gunIds))
                    Interactions = initialiseInteractionsStruct(Config);
                    warning(['Mismatch in the air gun identifiers between '...
                        '''Config'' and ''Interactions'' structures. '...
                        'Bubble interactions are neglected'])
                end
            else
                Interactions = initialiseInteractionsStruct(Config);
                warning(['Non-supported value for PROPERTY = '...
                    '''Interactions''. Bubble interactions are neglected'])
            end
            
        case 'calibration' % calibration parameters
            Calibration = values{m};
            isValid = false;
            if isstruct(Calibration)
                isValid = true; % true for valid calibration structure
                if ~isCalibrationStruct(Calibration)
                    isValid = false;
                    Calibration = initialiseCalibrationStruct(Config);
                end
            end
                
            if ~isValid
                warning(['Non-valid value for PROPERTY = ''Calibration'' '...
                    'The following values will be used for the fields in '...
                    'the ''Calibration'' structure: ETA = %0.3e, '...
                    'ALPHA = %0.3e, GAMMA = %0.3e'],Calibration.eta,...
                    Calibration.alpha,Calibration.gamma)
            end
            
        case 'excludephysics' % effective hydrostatic pressure
            options = values{m};
            validOptions = {'mu','zb'};
            if ~isempty(options)
                if ischar(options)
                    data = textscan(lower(options),'%s');
                    data = data{1};
                    isValid = ismember(data,validOptions);
                    excludePhysics = data(isValid);
                    if any(~isValid)
                        excludePhysics = '';
                        warning(['One or more character strings for the '...
                            'PROPERTY = ''ExcludePhysics'' are not supported'])
                    end
                else
                    warning(['The input for PROPERTY = ''ExcludePhysics'' must '...
                        'be a character vector'])
                end 
            end
            
        case 'displayprogress'
            displayProgress = values{m};
            if ~islogical(displayProgress) && ~any(displayProgress == [0 1])
                displayProgress = false;
                warning(['Non-supported value for PROPERTY = '...
                    '''DisplayProgress''. A value of 0 will be used'])
            end
    end
end

% Error Control (Dynamic Pressure in 'Interactions' Structure)
wrongDynamicPressure = false;
nGuns = length(gunIds); % number of selected air guns
for q = 1:nGuns
    ps = Interactions(q).dynamicPressure;
    if ~isnumeric(ps) || ~isvector(ps) || length(ps) ~= nPoints
        Interactions(q).dynamicPressure = zeros(1,nPoints);
        wrongDynamicPressure = wrongDynamicPressure || true;
    end
end
if wrongDynamicPressure
    warning(['The ''dynamicPressure'' field in the INTERACTIONS '...
        'structure must be a numeric vector of length '...
        'ROUND(SAMPLERATE * DURATION + 1). The hydrostatic pressure '...
        'will be used instead']) 
end
