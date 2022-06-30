function [targetAzimuth,targetElevation,spectrumType,metric,backprop] = ...
    plotSpectrum_ec(Radiation,varargin_cell)

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
validProperties = {'targetazimuth','targetelevation','spectrumtype',...
    'metric','backprop'};
properties = lower(varargin(1:2:end));
values = varargin(2:2:end);
if any(~ismember(properties,validProperties))
    error('One or more input properties are not recognised')
end

% Default Input Values
targetAzimuth = 0;
targetElevation = 90;
spectrumType = 'psd';
metric = 'rms';
backprop = true;

% Extract and Verify Input Values
nPairs = (nargin - nFixArg)/2; % number of (PROPERTY,VALUE) pairs
for m = 1:nPairs
    property = properties{m};
    switch property % populate with more properties if needed              
        case 'targetazimuth'
            targetAzimuth = values{m};
            if ~isempty(targetAzimuth) && (~isnumeric(targetAzimuth) ...
                    || ~isscalar(targetAzimuth) ...
                    || targetAzimuth < -180 || targetAzimuth > 180)
                targetAzimuth = 0;
                warning(['Non-valid value for PROPERTY = '...
                    '''TargetAzimuth''. A %0.0f deg target azimuth '...
                    'will be used'],targetAzimuth)
            end
        case 'targetelevation'
            targetElevation = values{m};
            if ~isempty(targetElevation) && (~isnumeric(targetElevation) ...
                    || ~isscalar(targetElevation) ...
                    || targetElevation < 0 || targetElevation > 90)
                targetElevation = 90;
                warning(['Non-valid value for PROPERTY = '...
                    '''TargetAzimuth''. A %0.0f deg target azimuth '...
                    'will be used'],targetElevation)
            end     
        case 'spectrumtype'
            spectrumType = values{m};
            if ~ischar(spectrumType) ...
                    || ~any(strcmpi(spectrumType,{'psd','cpb'}))
                spectrumType = 'psd';
                warning(['Non-valid value for PROPERTY = ''SpectrumType''. '...
                    'The %s spectrum will be plotted'],upper(spectrumType))
            end
        case 'metric'
            metric = values{m};
            if ~ischar(metric) || ~any(strcmpi(metric,{'rms','sel'}))
                metric = 'rms';
                warning(['Non-valid value for PROPERTY = ''Metric''. '...
                    'The %s spectrum will be plotted'],upper(metric))
            end
        case 'backprop'
            backprop = values{m};
            if isempty(backprop) || (~isnumeric(backprop) ...
                    && ~islogical(backprop)) || ~any(backprop == [0 1])
                backprop = true;
                warning(['Non-valid value for PROPERTY = ''Backprop''. '...
                    'The source level spectrum will be plotted'])
            end 
    end
end
