function Interactions = initialiseInteractionsStruct(varargin)

narginchk(0,1);

Interactions_empty = struct('gunId',[],'gunIds_sou',[],'gunIds_all',[],...
    'dynamicPressure',[]);

if nargin
    Config = varargin{1};
    if isConfigStruct(Config)
        % Load Parameters
        fs = Config.sampleRate;
        tmax = Config.duration;
        gunIds = Config.gunIds;
        gunIds_all = [Config.Array.gunId];
        
        % Populate Interactions Structure
        nPoints = round(fs*tmax + 1);
        nGuns = length(gunIds);
        for q = 1:nGuns
            Interactions(q) = Interactions_empty;
            Interactions(q).gunId = gunIds(q);
            Interactions(q).gunIds_sou = gunIds(gunIds ~= gunIds(q));
            Interactions(q).gunIds_all = gunIds_all;
            Interactions(q).dynamicPressure = zeros(1,nPoints);
        end
    else
        Interactions = Interactions_empty;
    end
else
    Interactions = Interactions_empty;
end