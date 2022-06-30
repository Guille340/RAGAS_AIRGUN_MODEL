% Bubble = polytropicDynamics(Config,varargin)
%
% DESCRIPTION
% Computes the motion parameters of the bubble produced by the seismic air
% gun array and environment conditions specified in the input CONFIG structure
% using the polytropic relation to establish a relationship between pressure
% and radius (Gilmore, 1952; Ziolkowsky, 1970). The optional input structures 
% INTERACTIONS and CALIBRATION contain information for the calculation of 
% bubble dynamics with interactions (Ziolkowsky et al., 1982) and calibration 
% parameters other than the default. The function does not support clustered 
% guns or GI guns.
%
% POLYTROPICDYNAMICS.m returns a BUBBLE structure containing the input
% configuration data (substructures 'Config', 'Array', 'Calibration', 
% 'Interactions') and the bubble motion parameters (substructure 'Dynamics')
% for each air gun in the array.
%
% The BUBBLE structure can be later used by function BUBBLERADIATION.m to
% generate the acoustic signature for the array at any specified receiver
% location.
%
% INPUT ARGUMENTS
% - Config: configuration structure generated from the JSON configuration
%   file stored in <ROOT.DATAFOLDER> directory (config.json). For details
%   about the content of this structure see INITIALISECONFIGSTRUCT.m.
%
% INPUT PROPERTIES
% - 'Interactions': bubble interactions structure generated with function
%   BUBBLEINTERACTIONS. The structure contains the effective acoustic pressure 
%   produced on each bubble by the surrounding nearby bubbles. The structure
%   also includes the identifiers of the various interacting guns. For details
%   about the content of this structure see INITIALISEINTERACTIONSSTRUCT.m.
% - 'Calibration': calibration structure generated with function
%   INITIALISECALIBRATIONSTRUCT.m. The structure contains the calibration
%   parameters for the bubble dynamics model ('eta','alpha','gamma') For 
%   details about its content see the aforementioned function.
% - 'ExcludePhysics': character vector used to specify what physical processes
%   are to be excluded from the simulations. The physical processes are 
%   referenced as two-character strings. More than one process can be included
%   in the character vector by adding one-space separation. These are the
%   physical processes available for exclusion in POLYTROPICDYNAMICS.m
%   # 'mu': effective viscosity. The effective viscosity accounts for the 
%      increase in the dynamic viscosity of the water due to a turbulent flow.
%   # 'zb': raise of air bubble towards the water surface. The raise speed
%     of large bubbles from air guns is relatively quick (~2 m/s) and its 
%     effect can be aparent after ~100 ms.
% - 'DisplayProgress': TRUE for displaying the progress of the simulations in
%    the command window.
%
% OUTPUT ARGUMENTS
% - Bubble: structure containing the input configuration information and 
%   predicted time-depended bubble motion parameters (e.g. pressure, velocity, 
%   radius, temperature and mass)
%
% FUNCTION CALL
% 1. Bubble = polytropicDynamics(Config)
% 2. Bubble = polytropicDynamics(...,PROPERTY,VALUE)
%    Properties: 'Interactions', 'Calibration', 'ExcludePhysics',
%    'DisplayProgress'
% 
% FUNCTION DEPENDENCIES
% - polytropicDynamics_ec
% - initialiseBubbleStruct
% - surfaceTensionWater
% - dynamicViscosityWater
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% REFERENCES
% Gilmore, F. R. (1952) "The growth or collapse of a spherical bubble in a 
% viscous compressible liquid", California Institute of Technology, report 26-4
%
% Ziolkowski, A. (1970). "A method for calculating the output pressure 
% waveform from an air gun”, Geophys. J. R. Astr. Soc. 21(2), 137-161
%
% Ziolkowski,, A., Parkes, G., Hatton, L. and Haugland, T. (1982). "The 
% signature of an air gun array: computation from near-field measurements 
% including interactions", Geophysics 47(10), 1413-1421
%
% See also MULTIPHYSICSDYNAMICS

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function Bubble = polytropicDynamics(Config,varargin)

narginchk(1,9)

% ERROR CONTROL
[Interactions,Calibration,excludePhysics,displayProgress] ...
    = polytropicDynamics_ec(Config,varargin);

% GENERAL PARAMETERS
tmax = Config.duration; % signature duration [s]
fs = Config.sampleRate; % sampling frequency for returned parameters [Hz]
fsp = round(1e4/fs)*fs; % sampling frequency for parameter processing [Hz]
rhow_inf = Config.waterDensity; % unperturbed density of sea water [kg m-3]
cw_inf = Config.waterSoundSpeed; % unperturbed speed of sound in sea water [m s-1]
Tw = Config.waterTemperature; % temperature of sea water [K]
Sw = Config.waterSalinity; % salinity of sea water [ppm]
gunIds = Config.gunIds; % identifiers of selected air guns
bubbleFormula = Config.bubbleFormula; % equation for bubble acceleration
interaction = ~all(~[Interactions.dynamicPressure]); % true for interactions

% UNIT CONVERSION FACTORS
atm2Pa = 101325; % atmosphere to Pascal
psi2Pa = 6894.75729; % pound per square inch to Pascal
cuin2cum = 0.000016387064; % cubic inch to cubic metre

% CONSTANTS
g = 9.80665; % acceleration of gravity [m s-2]
cpa = 1005; % specific heat capacity of air at constant pressure [J kg-1 K-1]
cva = 718; % specific heat capacity of air at constant volume [J kg-1 K-1]
Ra = cpa - cva; % gas constant for dry air [J kg-1 K-1]
n = 7; % water constant [-]
B = 3000*atm2Pa; % water constant [Pa]
patm = 1*atm2Pa; % atmospheric pressure [Pa]
Dt = 1/fsp; % time resolution step

% EQUATIONS OF BUBBLE MOTION
switch lower(bubbleFormula)
    case 'gilmore'
        bubbleAcceleration = @(rb,ub,pb,hb,dpb,dhb,cw,rhow,pw) ...
            ((1+ub/cw)*hb + rb/cw*(1-ub/cw)*dhb -3/2*ub^2*(1-ub/(3*cw))) ...
                / (rb*(1-ub/cw));

    case 'keller-kollodner'
        bubbleAcceleration = @(rb,ub,pb,hb,dpb,dhb,cw,rhow,pw) ...
            ((1+ub/cw)*((pb-pw)/rhow) + rb/(rhow*cw)*dpb ...
            -3/2*ub^2*(1-ub/(3*cw))) / (rb*(1-ub/cw));

    case 'herring'
        bubbleAcceleration = @(rb,ub,pb,hb,dpb,dhb,cw,rhow,pw) ...
            ((pb-pw)/rhow + rb/(rhow*cw)*dpb -3/2*ub^2*(1-4*ub/(3*cw))) ...
            / (rb*(1 - 2*ub/cw));

    case 'rayleigh-plesset'
        bubbleAcceleration = @(rb,ub,pb,hb,dpb,dhb,cw,rhow,pw) ...
            ((pb-pw)/rhow - 3/2*ub^2) / rb;    
end

% DISPLAY PROGRESS 
if displayProgress
    fprintf('\nPROCESSING BUBBLE DYNAMICS\n')
end

% PROCESSING
Bubble = initialiseBubbleStruct();
nGuns = length(gunIds);
for q = 1:nGuns
    
    % AIR GUN INDICES
    gunId = gunIds(q);
    iGun_con = find([Config.Array.gunId] == gunId);
    iGun_cal = find([Calibration.gunId] == gunId);
    iGun_int = find([Interactions.gunId] == gunId);  
    
    % DISPLAY PROGRESS
    if displayProgress
        fprintf('# Air gun ID = [%d] (%d/%d)\n',gunId,q,nGuns)
    end
    
    % EMPIRICAL CALIBRATION PARAMETERS
    eta = Calibration(q).eta; % air gun efficient
    alpha = Calibration(q).alpha; % attenuation exponent
    gamma = Calibration(q).gamma; % heat capacity ratio (cp/cv)

    % AIR GUN PARAMETERS
    pg = Config.Array(iGun_con).gunPressure; % pressure in airgun chambers (same for generator and injector) [psi]
    Vgg = Config.Array(iGun_con).gunGeneratorVolume; % volume of airgun chamber (generator) [in3] NOTE: for both G and GI guns
    z = Config.Array(iGun_con).gunZ; % airgun depth [m]
    ps = Interactions(iGun_int).dynamicPressure; % dynamic pressure from interactions

    % CONSTANTS AND CONVERSIONS
    ps = ps(:)'; % make it a row vector
    rg = 0; % airgun radius [m] (rg = 0 ignores gun presence)
    Vgun = 4*pi/3*rg^3; % effective volume of the gun (Vgun = 0 ignores gun pressence) [m3]
    pgt = pg*psi2Pa; % total pressure in the airgun chamber [Pa]
    Vgg = Vgg*cuin2cum; % volume in the airgun chamber (generator) [m3]

    % INITIAL CONDITIONS
    % Water Parameters
    pw_inf = patm + rhow_inf*g*z; % initial non-disturbed ambient pressure (= hydrostatic pressure)[Pa]
    pw0 = pw_inf + ps(1); % initial disturbed ambient pressure (affected by bubble raise and effective pressure) [Pa]
    rhow0 = rhow_inf; % initial disturbed density of water at bubble wall [kg m-3]
    cw0 = cw_inf; % initial disturbed sound speed in water at bubble wall [m s-1]
    sigmaw = surfaceTensionWater(Tw,Sw); % surface tension of sea water at 1 atm [N m-1]
    muw = dynamicViscosityWater(Tw,Sw); % dynamic viscosity of sea water at 1 atm [kg m-1 s-1]
    muwe = muw; % effective dynamic viscosity of water to account for turbulent flow [W m-1 K-1]

    % Gun Parameters
    Tgg = Tw*(1+pgt/139e6); % initial temperature in the gun chamber (generator) [ºC]
    mgtg = pgt*Vgg/(Ra*Tgg); % total mass in the gun chamber (generator) [kg]

    % Bubble Equilibrium Radius
    rbeq = (3*eta*mgtg*Ra*Tw/(pw_inf*4*pi))^(1/3);

    % Bubble Parameters
    ub0 = 0; % initial bubble wall velocity [m s-1]
    hb0 = 0; % initial specific enthalpy [m2 s-2]
    zb0 = z; % initial bubble depth [m]
    Vb0 = Vgg; % initial bubble volume [m3]
    rb0 = (3*(Vb0+Vgun)/(4*pi))^(1/3); % initial bubble radius [m]
    pb0 = pw_inf; % used instead of pb0 = pgt + patm from Ziolkowsky (1970)

    % Bubble Pressure
    pa0 = pb0 + 2*sigmaw/rb0;

    % Other Parameters
    rbcu = rb0^3; % third power of bubble radius (accumulated)
    dpb = 0;
    dhb = 0;
    dub = 0;

    % INITIALISE OUTPUT VECTOR
    nPoints = round(tmax * fsp + 1); % number of points in the airgun signature
    t = (0:nPoints - 1)/fsp; % time axis [s]

    % Water Parameters
    rhow = [rhow0, zeros(1,nPoints-1)];
    cw = [cw0, zeros(1,nPoints-1)];
    pw = [pw0, zeros(1,nPoints-1)];

    % Bubble Parameters
    pb = [pb0, zeros(1,nPoints-1)];
    ub = [ub0, zeros(1,nPoints-1)];
    rb = [rb0, zeros(1,nPoints-1)];
    hb = [hb0, zeros(1,nPoints-1)];
    Vb = [Vb0, zeros(1,nPoints-1)];
    zb = [zb0, zeros(1,nPoints-1)];

    % Bubble Pressure
    pa = [pa0, zeros(1,nPoints-1)];

    for m = 1:nPoints-1
        
        % SET OLD VALUES
        dub_old = dub;
        pb0_old = pb0;
        hb0_old = hb0;

        % TIME-DIFFERENTIAL QUATITIES
        % Bubble Depth
        if ~ismember('zb',excludePhysics)
            dzb = 2*g*rbcu*Dt/rb0.^3; % differential variation of bubble depth with time [m s-1] 
            rbcu = rbcu + rb0^3; % faster than dzb = 2*g*sum(rb(1:m).^3)*Dt/rb0.^3; 
        else
            dzb = 0;
        end

        % Radial Acceleration of the Bubble
        dub = bubbleAcceleration(rb0,ub0,pb0,hb0,dpb,dhb,cw0,rhow0,pw0);
        d2ub = (dub - dub_old)/Dt; % time derivative of the particle acceleration in the bubble wall [m s-3]

        % NEXT VALUES  
        % Bubble Parameters
        ub0 = ub0 + dub*Dt + d2ub*Dt^2/2; % next bubble wall velocity [m s-1] NOTE!: 2nd derivative slightly reduces the gap between peaks and valleys  
        rb0 = rb0 + ub0*Dt + dub*Dt^2/2; % next bubble radius [m] NOTE!: 2nd derivative slightly reduces the gap between peaks and valleys 
        Vb0 = 4/3*pi*rb0^3 - Vgun; % next bubble volume [m3]  
        zb0 = zb0 - dzb*Dt; % next bubble depth [m]

        % Bubble Pressure
        pa0 = pw0*(rbeq/rb0)^(3*gamma);
        pten = 2*sigmaw/rb0; % pressure due to tension at bubble wall [Pa]
        pvis = 4*muwe*ub0/rb0; % pressure due to viscosity of water at bubble wall [Pa]
        pb0 = pa0 - pten - pvis; % pressure at the external face of the bubble wall [Pa]
        hb0 = (pw0+B)/rhow0 * n/(n-1) * (((pb0+B)/(pw0+B))^((n-1)/n) -1);  % next specific enthalpy (= enthalpy / mass) [J kg-1]

        % Water Parameters
        pw0 = patm + rhow_inf*g*zb0 + ps(m+1); % next disturbed hydrostatic pressure [Pa] NOTE: using pw0 (time dep.) introduces a small error, since dhb is calculated assumming pw0 = cst
        rhow0 = rhow_inf*(((pb0+B)/(pw0+B))^(1/n)); % perturbed water density [Kg m-3]
        cw0 = cw_inf*((pb0+B)/(pw0+B))^((n-1)/(2*n)); % perturbed water sound speed [m s-1]

        % Effective Viscosity
        if ~ismember('mu',excludePhysics)
            Re = rhow0*abs(ub0)*2*rb0/muw; % Reynolds number
            muwe = muw*(1 + 0.02*Re); % effective dynamic viscosity of water [W m-1 K-1]
        else
            muwe = muw;
        end

        % Other Parameters
        dpb = (pb0 - pb0_old)/Dt;
        dhb = (hb0 - hb0_old)/Dt;

        % BUILD OUTPUT VECTOR (next value)
        % Water Parameters
        rhow(m+1) = rhow0;
        cw(m+1) = cw0;
        pw(m+1) = pw0; 

        % Bubble Parameters
        pb(m+1) = pb0;
        ub(m+1) = ub0;
        rb(m+1) = rb0;
        Tb(m+1) = NaN;
        mb(m+1) = NaN;
        hb(m+1) = hb0;
        Vb(m+1) = Vb0;
        zb(m+1) = zb0;

        % Bubble Pressure
        pf(m+1) = pa0;
        pa(m+1) = pa0;
        pv(m+1) = NaN;
        pag(m+1) = pa0;
        pai(m+1) = NaN;

        % Bubble Mass
        ma(m+1) = NaN;
        mv(m+1) = NaN;
        mag(m+1) = NaN;
        mai(m+1) = NaN;

        % Gun Parameters
        pgg(m+1) = NaN;
        pgi(m+1) = NaN;
        mgg(m+1) = NaN;
        mgi(m+1) = NaN;
        Tgg(m+1) = NaN;
        Tgi(m+1) = NaN;
    end

    % DAMPED PARAMETERS
    rb = exp(-alpha*t).*(rb - rbeq) + rbeq;
    ub = exp(-alpha*t).*(ub - alpha*(rb - rbeq));
    vb = 2*g*cumsum(rb.^3*Dt)./rb.^3; % bubble raise velocity [m/s]
    zb = z - cumsum(vb*Dt); % bubble depth [m]
    pw = patm + rhow*g.*zb + ps;
    Re = rhow.*abs(ub)*2.*rb/muw;
    muwe = muw*(1 + 0.02*Re);
    pten = 2*sigmaw./rb;
    pvis = 4*muwe.*ub./rb;
    pa =  pw.*(rbeq./rb).^(3*gamma);
    pb = pa - pten - pvis; % pressure at the external face of the bubble wall [Pa]
    rhow = rhow_inf*(((pb+B)./(pw+B)).^(1/n)); % perturbed water density [Kg m-3]
    cw = cw_inf*((pb+B)./(pw+B)).^((n-1)/(2*n)); % perturbed sound speed [m s-1]
    Vb = 4/3*pi*rb.^3 - Vgun; 
    hb = (pw+B)./rhow * n/(n-1) ...
        .* (((pb+B)./(pw+B)).^((n-1)/n) -1);  % next spec

    % SET VALUES TO NaN AFTER THE BUBBLE REACHES THE SURFACE
    ind = find(zb < 0,1,'first'); % index of first non-valid value

    % Water Parameters
    rhow(ind:nPoints) = NaN;
    cw(ind:nPoints) = NaN;
    pw(ind:nPoints) = NaN;

    % Bubble Parameters
    pb(ind:nPoints) = NaN;
    ub(ind:nPoints) = NaN;
    rb(ind:nPoints) = NaN;
    Tb(ind:nPoints) = NaN;
    mb(ind:nPoints) = NaN;
    hb(ind:nPoints) = NaN;
    Vb(ind:nPoints) = NaN;
    zb(ind:nPoints) = NaN;

    % Bubble Pressure
    pf(ind:nPoints) = NaN;
    pa(ind:nPoints) = NaN;
    pv(ind:nPoints) = NaN;
    pag(ind:nPoints) = NaN;
    pai(ind:nPoints) = NaN;

    % Bubble Mass
    ma(ind:nPoints) = NaN;
    mv(ind:nPoints) = NaN;
    mag(ind:nPoints) = NaN;
    mai(ind:nPoints) = NaN;

    % Gun Parameters
    pgg(ind:nPoints) = NaN;
    pgi(ind:nPoints) = NaN;
    mgg(ind:nPoints) = NaN;
    mgi(ind:nPoints) = NaN;
    Tgg(ind:nPoints) = NaN;
    Tgi(ind:nPoints) = NaN;

    % DOWNSAMPLE PARAMETERS TO TARGET SAMPLING RATE
    decimFactor = round(fsp/fs); % decimation factor
    t = t(1:decimFactor:end); % time axis [s]

    % Water Parameters
    rhow = rhow(1:decimFactor:end);
    cw = cw(1:decimFactor:end);
    pw = pw(1:decimFactor:end);

    % Bubble Parameters
    pb = pb(1:decimFactor:end);
    ub = ub(1:decimFactor:end);
    rb = rb(1:decimFactor:end);
    Tb = Tb(1:decimFactor:end);
    mb = mb(1:decimFactor:end);
    hb = hb(1:decimFactor:end);
    Vb = Vb(1:decimFactor:end);
    zb = zb(1:decimFactor:end);

    % Bubble Pressure
    pf = pf(1:decimFactor:end);
    pa = pa(1:decimFactor:end);
    pv = pv(1:decimFactor:end);
    pag = pag(1:decimFactor:end);
    pai = pai(1:decimFactor:end);

    % Bubble Mass
    ma = ma(1:decimFactor:end);
    mv = mv(1:decimFactor:end);
    mag = mag(1:decimFactor:end);
    mai = mai(1:decimFactor:end);

    % Gun Parameters
    pgg = pgg(1:decimFactor:end);
    pgi = pgi(1:decimFactor:end);
    mgg = mgg(1:decimFactor:end);
    mgi = mgi(1:decimFactor:end);
    Tgg = Tgg(1:decimFactor:end);
    Tgi = Tgi(1:decimFactor:end);

    % Populate BUBBLE Structure
    Bubble.Config = struct(...
        'gunIds',gunIds,...
        'duration',Config.duration,...
        'sampleRate',Config.sampleRate,...
        'waterDensity',Config.waterDensity,...
        'waterSoundSpeed',Config.waterSoundSpeed,...
        'waterTemperature',Config.waterTemperature,...
        'waterSalinity',Config.waterSalinity,...
        'bubbleMethod',Config.bubbleMethod,...
        'bubbleFormula',Config.bubbleFormula,...
        'interaction',Config.interaction);

    Bubble.Array(q) = struct(...
        'gunId',Config.Array(iGun_con).gunId,...
        'gunPressure',Config.Array(iGun_con).gunPressure,...
        'gunGeneratorVolume',Config.Array(iGun_con).gunGeneratorVolume,...
        'gunInjectorVolume',Config.Array(iGun_con).gunInjectorVolume,...
        'gunInjectorTime',Config.Array(iGun_con).gunInjectorTime,...
        'gunRadius',Config.Array(iGun_con).gunRadius,...
        'gunX',Config.Array(iGun_con).gunX,...
        'gunY',Config.Array(iGun_con).gunY,...
        'gunZ',Config.Array(iGun_con).gunZ);
    
    Bubble.Calibration(q) = Calibration(iGun_cal); %#ok<*FNDSB>

    Bubble.Interactions(q) = Interactions(iGun_int);
   
    Bubble.Dynamics(q) = struct('gunId',[],'interaction',interaction,...
        'timeAxis',t,'bubblePressure',pb,'bubbleVelocity',ub,...
        'bubbleRadius',rb,'bubbleEquilibriumRadius',rbeq,...
        'bubbleTemperature',Tb,'bubbleMass',mb,'bubbleEnthalpy',hb,...
        'bubbleVolume',Vb,'bubbleDepth',zb,'bubblePressureInternal',pf,...
        'bubblePressureAir',pa,'bubblePressureVapour',pv,...
        'bubblePressureAirGenerator',pag,'bubblePressureAirInjector',pai,...
        'bubbleMassAir',ma,'bubbleMassVapour',mv,...
        'bubbleMassAirGenerator',mag,'bubbleMassAirInjector',mai,...
        'gunPressureGenerator',pgg,'gunPressureInjector',pgi,...
        'gunMassAirGenerator',mgg,'gunMassAirInjector',mgi,...
        'gunTemperatureGenerator',Tgg,'gunTemperatureInjector',Tgi,...
        'waterDensity',rhow,'waterSoundSpeed',cw,'waterPressure',pw);
end
