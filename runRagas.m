function runRagam()

% Add Libraries
addpath(genpath(fullfile(pwd,'Core')))
addpath(genpath(fullfile(pwd,'Libraries')))

% Read Configuration and Plotting Data
Root = jsondecode(fileread(fullfile(pwd,'root.json'))); % root folder
Config = readConfig(Root); % verify and update configuration parameters

% Load Configuration Parameters
bubbleMethod = Config.bubbleMethod; % bubble motion modelling method

% Create Calibration Structure
Calibration = initialiseCalibrationStruct(Config);

% Process Bubble Dynamics (without Interactions)
Interactions = initialiseInteractionsStruct(Config);
if strcmp(bubbleMethod,'multiphysics')
    Bubble = multiphysicsDynamics(Config,'Interactions',Interactions,...
        'Calibration',Calibration,'ExcludePhysics','vg',...
        'DisplayProgress',true);
else % bubbleMethod = 'polytropic'
    Bubble = polytropicDynamics(Config,'Interactions',Interactions,...
        'Calibration',Calibration,'ExcludePhysics','',...
        'DisplayProgress',true);
end
filePath = fullfile(Root.dataFolder,'Bubble_Isolated.mat');
save(filePath,'-struct','Bubble')

% Process Bubble Dynamics (with Interactions)
Interactions = bubbleInteractions(Config,Bubble,'DisplayProgress',true);
if strcmp(bubbleMethod,'multiphysics')
    Bubble_int = multiphysicsDynamics(Config,'Interactions',Interactions,...
        'Calibration',Calibration,'ExcludePhysics','vg',...
        'DisplayProgress',true);
else % bubbleMethod = 'polytropic'
    Bubble_int = polytropicDynamics(Config,'Interactions',Interactions,...
        'Calibration',Calibration,'DisplayProgress',true);
end
filePath = fullfile(Root.dataFolder,'Bubble_Interacting.mat');
save(filePath,'-struct','Bubble_int')
