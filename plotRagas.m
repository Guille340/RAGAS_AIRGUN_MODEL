
function plotRagam()

% Add Libraries
addpath(genpath(fullfile(pwd,'Core')))
addpath(genpath(fullfile(pwd,'Libraries')))

% Read Configuration and Plotting Data
Root = jsondecode(fileread(fullfile(pwd,'root.json'))); % root folder
Config = readConfig(Root); % verify and update configuration parameters
Plot = readPlot(Root);

% Load Configuration Parameters
gunIds = Config.gunIds; % air gun identifiers to be processed
interaction = Plot.interaction;
ghost = Plot.radiationGhost; % ghost reflection mode
radius = Plot.souToRecDistance;
frequencies = Plot.frequency;
targetAzim = Plot.souToRecAzimuth;
targetElev = Plot.souToRecElevation;
   
% Vertical Directivity
nElev = 18;
azim = [targetAzim, targetAzim - 180];
elev = 0:90/(nElev - 1):90;
[x_arr,y_arr,z_arr] = arrayCentre(Config);
[x_rec,y_rec,z_rec] = gridSphere(radius,[x_arr,y_arr,z_arr],azim,elev);
if interaction
    filePath = fullfile(Root.dataFolder,'Bubble_Interacting.mat');
    Bubble = load(filePath);
    Radiation = bubbleRadiation(Bubble,[x_rec y_rec z_rec],...
        'GunIdentifiers',gunIds,'Ghost',ghost,'DisplayProgress',true);
else 
    filePath = fullfile(Root.dataFolder,'Bubble_Isolated.mat');
    Bubble = load(filePath);
    Radiation = bubbleRadiation(Bubble,[x_rec y_rec z_rec],...
        'GunIdentifiers',gunIds,'Ghost',ghost,'DisplayProgress',true);
end
h_ver = plotVerticalDirectivity(Radiation,'frequencies',frequencies,...
    'targetazimuth',targetAzim);

% Horizontal Directivity
nAzim = 72;
azim = 0:360/(nAzim - 1):360;
elev = targetElev;
[x_arr,y_arr,z_arr] = arrayCentre(Bubble);
[x_rec,y_rec,z_rec] = gridSphere(radius,[x_arr,y_arr,z_arr],azim,elev);
if interaction
    filePath = fullfile(Root.dataFolder,'Bubble_Interacting.mat');
    Bubble = load(filePath);
    Radiation = bubbleRadiation(Bubble,[x_rec y_rec z_rec],...
        'GunIdentifiers',gunIds,'Ghost',ghost,'DisplayProgress',true);
else 
    filePath = fullfile(Root.dataFolder,'Bubble_Isolated.mat');
    Bubble = load(filePath);
    Radiation = bubbleRadiation(Bubble,[x_rec y_rec z_rec],...
        'GunIdentifiers',gunIds,'Ghost',ghost,'DisplayProgress',true);
end
h_hor = plotHorizontalDirectivity(Radiation,'frequencies',frequencies,...
    'targetelevation',targetElev);

% Source Level Spectrum
h_psd = plotSpectrum(Radiation,'TargetAzimuth',targetAzim,...
    'TargetElevation',targetElev,'SpectrumType','psd','Metric','rms');
h_cpb = plotSpectrum(Radiation,'TargetAzimuth',targetAzim,...
    'TargetElevation',targetElev,'SpectrumType','cpb','Metric','rms');

% Array Signature
h_sig = plotSignature(Radiation,'TargetAzimuth',targetAzim,...
    'TargetElevation',targetElev,'backprop',true);

% Save Figures
figDir = fullfile(Root.dataFolder,'Figures');
mkdir(figDir)
figName_ver = sprintf('Vertical Directivity (azim = %0.0f)',targetAzim);
savefig(h_ver,fullfile(figDir,figName_ver))
print(h_ver,fullfile(figDir,figName_ver),'-dpng','-r200')
fprintf('\n# Saving vertical directivity ...\n')

figName_hor = sprintf('Horizontal Directivity (elev = %0.0f)',targetElev);
savefig(h_hor,fullfile(figDir,figName_hor))
print(h_hor,fullfile(figDir,figName_hor),'-dpng','-r200')
fprintf('# Saving horizontal directivity ...\n')

figName_psd = sprintf('PSD (azim = %0.0f, elev = %0.0f)',targetAzim,targetElev);
savefig(h_psd,fullfile(figDir,figName_psd))
print(h_psd,fullfile(figDir,figName_psd),'-dpng','-r200')
fprintf('# Saving power spectral density spectrum ...\n')

figName_cpb = sprintf('CPB (azim = %0.0f, elev = %0.0f)',targetAzim,targetElev);
savefig(h_cpb,fullfile(figDir,figName_cpb))
print(h_cpb,fullfile(figDir,figName_cpb),'-dpng','-r200')
fprintf('# Saving band spectrum ...\n')

figName_sig = sprintf('Array Signature (azim = %0.0f, elev = %0.0f)',...
    targetAzim,targetElev);
savefig(h_sig,fullfile(figDir,figName_sig))
print(h_sig,fullfile(figDir,figName_sig),'-dpng','-r200')
fprintf('# Saving array signature ...\n')


