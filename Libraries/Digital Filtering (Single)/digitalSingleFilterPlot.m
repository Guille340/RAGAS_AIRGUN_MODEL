%  hfig = DIGITALSINGLEFILTERPLOT(DigitalFilter)
%
%  DESCRIPTION
%  Plots the decimators (red, if applicable) and filter (magenta) defined in 
%  the input structure DIGITALFILTER (see DIGITALSINGLEFILTERDESIGN). The 
%  markers on the decimator and digital filter curves designate the half-power 
%  (-3 dB) frequencies.
% 
%  INPUT VARIABLES
%  - DigitalFilter: structure containing the filtering object and information
%    of a filter generated with DIGITALFILTERDESIGN.
%   
%  OUTPUT VARIABLES
%  - hfig: figure handle
%
%  FUNCTION CALL
%  hfig = digitalSingleFilterPlot(DigitalFilter)
%
%  FUNCTION DEPENDENCIES
%  - None
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  EXAMPLE
%  % 1) Configuration Data
%  fs = 44100;
%  bpo = 3;
%  freqCutoff = [20 1e3];
%  filtOrder = 10;
%
%  % 2) Filter Bank Design
%  DigitalFilter = digitalSingleFilterDesign(fs,freqCutoff,'FilterOrder',...
%      filtOrder)
%
%  % 3) Plot Filter Bank
%  hfig = digitalSingleFilterPlot(DigitalFilter)
%
%  See also DIGITALFILTERDESIGN, DIGITALFILTER

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  06 Aug 2021

function hfig = digitalSingleFilterPlot(DigitalFilter)
   
% Error Control (check filter structure)
genfieldsValid = {'sampleRate','freqCutoff','decimator','target'};
decfieldsValid = {'filterType','filterResponse','filterOrder',...
    'filterObject','halfPowerFreqn','decimationFactor','sampleRate',...
    'halfPowerFreq'};
bpafieldsValid = {'filterType','filterResponse','filterOrder',...
    'filterObject','halfPowerFreqn1','halfPowerFreqn2',...
    'sampleRate','halfPowerFreq1','halfPowerFreq2'};

genfields = fieldnames(DigitalFilter)';
if isequal(genfieldsValid,genfields)
    decfields = fieldnames(DigitalFilter.decimator)';
    bpafields = fieldnames(DigitalFilter.target)';
    if ~isequal(decfieldsValid,decfields) || ~isequal(bpafieldsValid,bpafields)
        error('The input argument DIGITALFILTER is not a valid filter structure')
    end
else
    error('The input argument DIGITALFILTER is not a valid filter structure')
end

% Load Parameters
decimator = DigitalFilter.decimator;
target = DigitalFilter.target;
freqLims = DigitalFilter.freqCutoff;
fs = DigitalFilter.sampleRate;
filtOrder = target.filterOrder;
decimFilter = decimator.filterObject;
decimFactor = decimator.decimationFactor;
fdn1 = decimator.halfPowerFreqn;
fbn1 = target.halfPowerFreqn1;
fbn2 = target.halfPowerFreqn2;

% Plot Filter Bank with ANSI Masks
hfig = [];
nDecimSteps = log2(decimFactor); % number of parent octave bands
if ~isempty(nDecimSteps) && nDecimSteps > 0
    hfig = figure;
    hold on
    nPointsInOctave = 50; % plotting points per fractional octave band
    fmin = freqLims(1)*2^-3; % lower frequency limit for calculations
    fmax = fs/2; % upper frequency limit for calculations
    nPoints = round(nPointsInOctave * log(fmax/fmin)/log(2)); % total number of plotting points for each curve
    fAxis = fmin * (fmax/fmin).^((0:nPoints-1)/(nPoints-1)); % frequency axis
    nDecimators = sum(~isnan(decimator.halfPowerFreq)); % number of decimation filters
    decimHandleIndex = nan(1,nDecimators);
    markHandleIndex = nan(1,nDecimators);
    plotCounter = 0;
    iDecim = 0;
    for m = 1:nDecimSteps
        % Plot Decimation Filter for Current Octave Band
        D = decimFactor/2^(nDecimSteps-m+1);    
        if D >= 2
            fr = 2*fs/D; % sampling rate (before downsampling by 2)
            fAxis = fAxis(fAxis <= fr/2); % limit frequency axis to Nyquist frequency
            decimResp = 20*log10(abs(freqz(decimFilter,fAxis*2/fr,2))); % frequency response of decimator

            hDecim = plot(fAxis,decimResp,'r','LineWidth',1); % plot decimation filter before decimation
            [~,iDecimMark] = min(abs(fAxis*2/fr - fdn1)); % axis index of half-power frequency of decimator
            plot(fdn1*fr/2,decimResp(iDecimMark),'r','LineWidth',1,...
                'Marker','diamond','MarkerFaceColor','w',...
                'MarkerEdgeColor','r','MarkerSize',5); % plot mark 2 (half-power frequency of decimator)
            plotCounter = plotCounter + 2; 
            iDecim = iDecim + 1;
            decimHandleIndex(iDecim) = plotCounter-2; % handle indices for decimator curve
            markHandleIndex(iDecim) = plotCounter-1; % handle indices for mark 1 (central freq of octave)
        end

        plotCounter = plotCounter + 1;
    end
    
    % Plot Target Filter for Current Octave Band
    fr = fs/D; % sampling rate (after downsampling by 2)
    fAxis = fAxis(fAxis <= fr/2); % limit frequency axis to Nyquist frequency
    targetResp = 20*log10(abs(freqz(target.filterObject,fAxis*2/fr,2))); % current bandpass filter response
    hTarget = plot(fAxis,targetResp,'m','LineWidth',1.5); % plot bandpass filter response
    [~,iTargetMark1] = min(abs(fAxis*2/fr - fbn1)); % axis index of half-power frequency of decimator
    [~,iTargetMark2] = min(abs(fAxis*2/fr - fbn2)); % axis index of half-power frequency of decimator
    plot(fbn1*fr/2,targetResp(iTargetMark1),'m','LineWidth',1,...
        'Marker','o','MarkerFaceColor','w','MarkerEdgeColor','m',...
        'MarkerSize',5); % plot half-power frequency 1 of target filter
    plot(fbn2*fr/2,targetResp(iTargetMark2),'m','LineWidth',1,...
        'Marker','o','MarkerFaceColor','w','MarkerEdgeColor','m',...
        'MarkerSize',5); % plot half-power frequency 2 of target filter

    % Plot Settings
    set(gca,'XScale','log')
    title({'Response of Single Digital Filter';...
        sprintf(['\\rm\\fontsize{10}IIR Order %d, f1 = %0.0f Hz, '...
        'f2 = %0.0f Hz'],filtOrder,fbn1*fr/2,fbn2*fr/2)})
    warning('off','MATLAB:legend:IgnoringExtraEntries')
    legend([hTarget hDecim],{'Digital Filter','Decimator'},'Location',...
        'NorthEast')
    axis([fmin fs/2 -80 10])
    xlabel('Frequency [Hz]')
    ylabel('Filter Gain [dB]')
    box on
    grid on
    set(gcf,'Units','Normalized','OuterPosition',[0.05 0.07 0.9 0.90])
else
    warning('Input DIGITALFILTER is an empty filter structure')
end
