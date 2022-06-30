function [gunIds,ghost,displayProgress] = bubbleRadiation_ec(Bubble,...
    receiverPositions,varargin_cell)

% Verify number of Input Arguments
nFixArg = 2;
varargin = varargin_cell;
nVarargin = length(varargin);
nargin = nFixArg + nVarargin;
if rem(nargin - nFixArg,2)
    error('Variable input arguments must come in pairs (PROPERTY,VALUE)')
end

% Error Control ('Bubble')
if ~isBubbleStruct(Bubble)
    error('Input BUBBLE is not a valid ''Bubble'' structure')
end

% Error Control ('receiverPositions')
if ~isnumeric(receiverPositions) || ~ismatrix(receiverPositions)
    error('Input RECEIVERPOSITIONS must be a numeric vector or matrix')
end

% Extract and Verify Input Properties
validProperties = {'gunidentifiers','ghost','displayprogress'};
properties = lower(varargin(1:2:end));
if any(~ismember(properties,validProperties))
    error('One or more input properties are not recognised')
end

% General Parameters
gunIds_all = Bubble.Config.gunIds;

% Default Input Values
gunIds = gunIds_all;
displayProgress = false;

% Extract and Verify Input Values
values = varargin(2:2:end);
nPairs = (nargin - nFixArg)/2; % number of (PROPERTY,VALUE) pairs
for m = 1:nPairs
    property = properties{m};
    switch property % populate with more properties if needed   
        case 'gunidentifiers'
            gunIds = values{m};
            if ~isnumeric(gunIds) || ~isvector(gunIds) ...
                    || any(~ismember(gunIds,gunIds_all))
                gunIds = gunIds_all;
                warning(['One or more values in PROPERTY = '...
                    '''GunIdentifiers'' are not valid air gun '...
                    'identifiers. Check the BUBBLE input'])
            end
            
        case 'ghost'
            ghost = values{m};
            if ~ismember(ghost,[0 1])
                ghost = true;
                warning(['Non-valid value for PROPERTY = ''Ghost''. '...
                    'The ghost reflection will be included (GHOST '...
                    '= %d'],ghost)
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
