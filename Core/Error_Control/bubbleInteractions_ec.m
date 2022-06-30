function [displayProgress] = bubbleInteractions_ec(Config,Bubble,varargin_cell)

% Verify number of Input Arguments
nFixArg = 2;
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

% Error Control ('Bubble')
if ~isBubbleStruct(Bubble)
    error('Input BUBBLE is not a valid ''Bubble'' structure')
end

% ERROR CONTROL ('Config' and 'Bubble' Matching)
if any(~ismember(Config.gunIds,[Bubble.Array.gunId])) ...
        || any(~ismember([Bubble.Array.gunId],Config.gunIds))
    error(['Mismatch in the air gun identifiers between ''Config'' and '...
        '''Bubble'' structures'])
end

% Extract and Verify Input Properties
validProperties = {'displayprogress'};
properties = lower(varargin(1:2:end));
if any(~ismember(properties,validProperties))
    error('One or more input properties are not recognised')
end

% Default Input Values
displayProgress = false;

% Extract and Verify Input Values
values = varargin(2:2:end);
nPairs = (nargin - nFixArg)/2; % number of (PROPERTY,VALUE) pairs
for m = 1:nPairs
    property = properties{m};
    switch property % populate with more properties if needed            
        case 'displayprogress'
            displayProgress = values{m};
            if ~islogical(displayProgress) && ~any(displayProgress == [0 1])
                displayProgress = false;
                warning(['Non-supported value for PROPERTY = '...
                    '''DisplayProgress''. A value of 0 will be used'])
            end
    end
end
