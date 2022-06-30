function [targetElevation,frequencies] = ...
    plotHorizontalDirectivity_ec(Radiation,varargin_cell)

% Verify number of Input Arguments
nFixArg = 1;
varargin = varargin_cell;
nVarargin = length(varargin);
nargin = nFixArg + nVarargin;
if rem(nargin - nFixArg,2)
    error('Variable input arguments must come in pairs (PROPERTY,VALUE)')
end

% Error Control ('Radiation')
if ~isRadiationStruct(Radiation)
    error('Input RADIATION is not a valid radiation structure')
end

% Extract and Verify Input Properties
validProperties = {'targetelevation','frequencies'};
properties = lower(varargin(1:2:end));
values = varargin(2:2:end);
if any(~ismember(properties,validProperties))
    error('One or more input properties are not recognised')
end

% Default Input Values
targetElevation = 45;
frequencies = [];

% Extract and Verify Input Values
nPairs = (nargin - nFixArg)/2; % number of (PROPERTY,VALUE) pairs
for m = 1:nPairs
    property = properties{m};
    switch property % populate with more properties if needed              
        case 'targetelevation'
            targetElevation = values{m};
            if ~isempty(targetElevation) && (~isnumeric(targetElevation) ...
                    || ~isscalar(targetElevation) ...
                    || targetElevation < 0 || targetElevation > 90)
                targetElevation = 0;
                warning(['Non-valid value for PROPERTY = '...
                    '''TargetElevation''. A %0.0f deg target elevation '...
                    'will be used'],targetElevation)
            end
            
        case 'frequencies'
            frequencies = values{m};
            if ~isempty(frequencies) && (~isnumeric(frequencies) ...
                    || ~isvector(frequencies))
                frequencies = [];
                warning(['Non-valid value for PROPERTY = ''Frequencies''. '...
                    'The broadband pattern will be plotted'])
            end
    end
end
