function Calibration = initialiseCalibrationStruct(varargin)

% option: bubbleMethod or Config structure.

narginchk(0,1);

Calibration_empty = struct('gunId',[],'taug',[],'taui',[],'beta',[],...
    'epsil',[],'kvis',[],'eta',[],'alpha',[],'gamma',[]);

if nargin
    Config = varargin{1};
    if isConfigStruct(Config)
        gunIds = Config.gunIds;
        nGuns = length(gunIds);
        if strcmpi(Config.bubbleMethod,'multiphysics')
            for q = 1:nGuns
                Calibration(q) = struct('gunId',gunIds(q),'taug',0.09115,...
                    'taui',0.09115,'beta',0.52,'epsil',22230,...
                    'kvis',0.02,'eta',0.8317,'alpha',[],'gamma',[]);
            end
        else
            for q = 1:nGuns
                Calibration(q) = struct('gunId',gunIds(q),'taug',[],...
                    'taui',[],'beta',[],'epsil',[],'kvis',[],...
                    'eta',0.8317,'alpha',10,'gamma',1.13);
            end
        end
    else
        Calibration = Calibration_empty;
    end
else
    Calibration = Calibration_empty;
end
