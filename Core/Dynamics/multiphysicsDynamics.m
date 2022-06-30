% Bubble = multiphysicsDynamics(Config,varargin)
%
% DESCRIPTION
% Computes the motion parameters of the bubble produced by the seismic air
% gun array and environment conditions specified in the input CONFIG structure
% using a multiphysics approach (Ziolkowsky, 1970; Laws et al, 1990). The 
% optional input structures INTERACTIONS and CALIBRATION contain information 
% for the calculation of bubble dynamics with interactions (Ziolkowsky et al.,
% 1982) and calibration parameters other than the default. The function 
% supports clustered guns (Laws et al, 1990) and generator-injector (GI) guns 
% (Landro, 1992). 
%
% MULTIPHYSICSDYNAMICS.m returns a BUBBLE structure containing the input
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
%   about the content of this structure see INITIALISEINTERACTIONSSTRUCT.m
% - 'Calibration': calibration structure generated with function
%   INITIALISECALIBRATIONSTRUCT.m. The structure contains the calibration
%   parameters for the bubble dynamics model ('taug','taui','beta','epsil',
%   'eta'). For details about its content see the aforementioned function.
% - 'ExcludePhysics': character vector used to specify what physical processes
%   are to be excluded from the simulations. The physical processes are 
%   referenced as two-character strings. More than one process can be included
%   in the character vector by adding one-space separation. These are the
%   physical processes available for exclusion in MULTIPHYSICSDYNAMICS.m
%   # 'jt': Joule-Thomson effect. This accounts for the temperature difference
%      between the air outside and inside the air gun during the throttling
%      process.
%   # 'mv': transfer rate of water vapour mass during condensation-vaporisation.
%   # 'qw': heat flow between bubble and water through conduction.
%   # 'qv': heat loss in the phase change from water to vapour. The energy
%      is proportional to the transfer rate of vapour mass and the latent 
%      heat of water.
%   # 'vg': air gun presence. The volume of the air gun affects the maximum
%      radius of the bubble. Including the air gun volume makes the code 
%      unstable. It is strongly recommended to exclude the air gun presence 
%      for the simulations.
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
% 1. Bubble = multiphysicsDynamics(Config)
% 2. Bubble = multiphysicsDynamics(...,PROPERTY,VALUE)
%    Properties: 'Interactions', 'Calibration', 'ExcludePhysics',
%    'DisplayProgress'
% 
% FUNCTION DEPENDENCIES
% - multiphysicsDynamics_ec
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
% Laws, R. M., Hatton, L. and Haartsen, M. (1990). "Computer modelling of
% clustered airguns", First Break 8(9), pp. 331-337
%
% Ziolkowski,, A., Parkes, G., Hatton, L. and Haugland, T. (1982). "The 
% signature of an air gun array: computation from near-field measurements 
% including interactions", Geophysics 47(10), 1413-1421
%
% Landrø, M. (1992). "Modelling of GI gun signatures", Geophysical Prospecting 
% 40(7), pp. 721-747
%
% See also POLYTROPICDYNAMICS

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function Bubble = multiphysicsDynamics(Config,varargin)

narginchk(1,9)

% ERROR CONTROL
[Interactions,Calibration,excludePhysics,displayProgress] ...
    = multiphysicsDynamics_ec(Config,varargin);

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

% INTERACTION PARAMETER
ps_all = [Interactions.dynamicPressure];
interaction = ~all(~ps_all(~isnan(ps_all))); % true for zero dynamic pressure

% UNIT CONVERSION FACTORS
atm2Pa = 101325; % atmosphere to Pascal
psi2Pa = 6894.75729; % pound per square inch to Pascal
cuin2cum = 0.000016387064; % cubic inch to cubic metre

% CONSTANTS
krou = 10; % roughness factor to account for the higher exposed area of a not perfectly spherical bubble (krou = 10 recommended)
kcor = 1; % correction factor for the equation of mass transfer of water vapour (kcor = 1 recommended)
g = 9.80665; % acceleration of gravity [m s-2]
cpa = 1005; % specific heat capacity of air at constant pressure [J kg-1 K-1]
cva = 718; % specific heat capacity of air at constant volume [J kg-1 K-1]
cpv = 4220; % specific heat capacity of water vapour at constant pressure (100 ºC) [J kg-1 K-1]
cvv = 3758.48; % specific heat capacity of water vapour at constant volume (100 ºC)[J kg-1 K-1]
Ra = cpa - cva; % gas constant for dry air [J kg-1 K-1]
Rv = cpv - cvv; % gas constant for water vapour [J kg-1 K-1]
n = 7; % water constant [-]
B = 3000*atm2Pa; % water constant [Pa]
patm = 1*atm2Pa; % atmospheric pressure [Pa]
kec = 0.04; % accomodation factor of evaporation-condensation
Lv = 2.26*1e6; % latent heat of water evaporation [J kg-1]
pc = 22.13*1e6; % pressure at the critical (saturated) state [Pa]
Tc = 647.31; % temperature at the critical (saturated) state [K]
kte = 6*1e-3; % coefficient of thermal expansion for dry air (20 ºC) [K-1]
pveq = pc*exp(Lv/(Rv*Tc)*(1-Tc/Tw)); % equilibrium vapour pressure [Pa]
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
    taug = Calibration(q).taug;
    taui = Calibration(q).taui;
    beta = Calibration(q).beta;
    eta = Calibration(q).eta;
    epsil = Calibration(q).epsil;
    kvis = Calibration(q).kvis;

    % AIR GUN PARAMETERS
    pg = Config.Array(iGun_con).gunPressure; % pressure in airgun chambers (same for generator and injector) [psi]
    Vgg = Config.Array(iGun_con).gunGeneratorVolume; % volume of airgun chamber (generator) [in3] NOTE: for both G and GI guns
    Vgi = Config.Array(iGun_con).gunInjectorVolume; % volume of airgun chamber (injector) [in3] NOTE: only for GI guns
    ti = Config.Array(iGun_con).gunInjectorTime * 1e-3; % start time for the injection [s] NOTE: only for GI guns
    rg = Config.Array(iGun_con).gunRadius; % airgun radius [m] NOTE: typically 14-14.5 cm, 10 cm for mini-guns. rg = 0 ignores gun presence
    z = Config.Array(iGun_con).gunZ; % airgun depth [m]
    ps = Interactions(iGun_int).dynamicPressure; % dynamic pressure from interactions

    % IGNORE GUN PRESENCE
    if ismember('vg',excludePhysics)
        rg = 0;
    end

    % CONSTANTS AND CONVERSIONS
    Vgun = 4*pi/3*rg^3; % effective volume of the gun (Vgun = 0 ignores gun pressence) [m3]
    pgt = pg*psi2Pa; % total pressure in the airgun chamber [Pa]
    Vgg = Vgg*cuin2cum; % volume in the airgun chamber (generator) [m3]
    Vgi = Vgi*cuin2cum; % volume in the airgun chamber (injector) [m3]

    % INITIAL CONDITIONS
    % Water Parameters
    pw_inf = patm + rhow_inf*g*z; % initial non-disturbed ambient pressure (= hydrostatic pressure)[Pa]
    pw0 = pw_inf; % initial non-disturbed ambient pressure (affected by bubble raise) [Pa]
    pe0 = pw_inf + ps(1); % initial disturbed ambient pressure (effective pressure)[Pa]
    rhow0 = rhow_inf; % initial disturbed density of water at bubble wall [kg m-3]
    cw0 = cw_inf; % initial disturbed sound speed in water at bubble wall [m s-1]
    sigmaw = surfaceTensionWater(Tw,Sw); % surface tension of sea water at 1 atm [N m-1]
    muw = dynamicViscosityWater(Tw,Sw); % dynamic viscosity of sea water at 1 atm [kg m-1 s-1]
    muwe = muw; % effective dynamic viscosity of water to account for turbulent flow [W m-1 K-1]

    % Bubble Parameters
    pb0 = pe0; % initial pressure in the external side of the bubble wall [Pa] 
    ub0 = 0; % initial bubble wall velocity [m s-1]
    hb0 = 0; % initial specific enthalpy [m2 s-2]
    zb0 = z; % initial bubble depth [m]
    Vb0 = Vgg; % initial bubble volume [m3]
    rb0 = (3*(Vb0+Vgun)/(4*pi))^(1/3); % initial bubble radius [m]
    Tb0 = Tw*(1+pgt/139e6); % initial temperature in the bubble [ºC] NOTE: used Tg0 instead of Tw!

    % Bubble Pressure
    pf0 = pb0 + 2*sigmaw/rb0; % initial pressure inside the bubble (pf0 = pa0 + pv0) [Pa]
    pag0 = pf0; % initial partial pressure of air in the bubble (from generator) [Pa]
    pai0 = 0; % initial partial pressure of air in the bubble (from injector) [Pa]
    pa0 = pag0 + pai0; % % initial partial pressure of air in the bubble [Pa]
    pv0 = 0; % initial partial pressure of water vapour in the bubble [Pa]

    % Bubble Mass
    mag0 = pag0*Vb0/(Ra*Tb0);  % initial mass of non-condensable gas (dry air) in the bubble (generator) [kg]
    mai0 = pai0*Vb0/(Ra*Tb0);  % initial mass of non-condensable gas (dry air) in the bubble (injector) [kg]
    ma0 = pa0*Vb0/(Ra*Tb0);  % initial mass of non-condensable gas (dry air) in the bubble (ma0 = mag0 + mai0) [kg]
    mv0 = 0; % initial mass of condensable gas (water vapour) in the bubble [kg]
    mb0 = ma0 + mv0; % initial mass of the gas mixture (air + vapour) in the bubble [kg]

    % Gun Parameters
    Tgg0 = Tw*(1+pgt/139e6); % initial temperature in the gun chamber (generator) [ºC]
    Tgi0 = Tw*(1+pgt/139e6); % initial temperature in the gun chamber (injector) [ºC]
    mgtg = pgt*Vgg/(Ra*Tgg0); % total mass in the gun chamber (generator) [kg]
    mgti = pgt*Vgi/(Ra*Tgi0); % total mass in the gun chamber (injector) [kg]
    mgg0 = mgtg - mag0; % initial mass in the gun chamber (generator) [kg] NOTE: used mgtg-mag0 instead of mgtg!
    mgi0 = mgti - mai0; % initial mass in the gun chamber (injector) [kg] NOTE: used mgti-mai0 instead of mgti!
    pgg0 = mgg0*Ra*Tgg0/Vgg; % initial pressure in the gun chamber (generator) [Pa] NOTE: used mgg0 instead of mgtg!
    pgi0 = mgi0*Ra*Tgi0/Vgi; % initial pressure in the gun chamber (injector) [Pa] NOTE: used mgi0 instead of mgti!

    % Set Partial Gun Mass and Volumes to EPS to avoid NaN Results
    pgg0 = max(pgg0,eps);
    pgi0 = max(pgi0,eps);
    mgg0 = max(mgg0,eps);
    mgi0 = max(mgi0,eps);

    % Other Parameters
    kjtg = Vgg/(mgg0*cpa)*(kte*Tgg0 - 1); % Joule-Thomson coefficient (generator) [K Pa-1]
    kjti = Vgi/(mgi0*cpa)*(kte*Tgi0 - 1); % Joule-Thomson coefficient (injector) [K Pa-1]
    rbcu = rb0^3; % third power of bubble radius (accumulated)
    dub = 0; % initial bubble wall acceleration [m s-2] 
    dmag = 0;
    dmai = 0;
    drhow = 0;
    dpw = 0;

    % INITIALISE OUTPUT VECTOR
    nPoints = round(tmax * fsp + 1); % number of points in the airgun signature
    t = (0:nPoints - 1)/fsp; % time axis [s]

    % Bubble Equilibrium Radius
    rbeq = (3*eta*(mgtg + mgti)*Ra*Tw/(pw0*4*pi))^(1/3);

    % Water Parameters
    rhow = [rhow0, zeros(1,nPoints-1)];
    cw = [cw0, zeros(1,nPoints-1)];
    pw = [pw0, zeros(1,nPoints-1)];
    pe = [pe0, zeros(1,nPoints-1)];

    % Bubble Parameters
    pb = [pb0, zeros(1,nPoints-1)];
    ub = [ub0, zeros(1,nPoints-1)];
    rb = [rb0, zeros(1,nPoints-1)];
    Tb = [Tb0, zeros(1,nPoints-1)];
    mb = [mb0, zeros(1,nPoints-1)];
    hb = [hb0, zeros(1,nPoints-1)];
    Vb = [Vb0, zeros(1,nPoints-1)];
    zb = [zb0, zeros(1,nPoints-1)];

    % Bubble Pressure
    pf = [pf0, zeros(1,nPoints-1)];
    pa = [pa0, zeros(1,nPoints-1)];
    pag = [pag0, zeros(1,nPoints-1)];
    pai = [pai0, zeros(1,nPoints-1)];
    pv = [pv0, zeros(1,nPoints-1)];

    % Bubble Mass
    ma = [ma0, zeros(1,nPoints-1)];
    mag = [mag0, zeros(1,nPoints-1)];
    mai = [mai0, zeros(1,nPoints-1)];
    mv = [mv0, zeros(1,nPoints-1)];

    % Gun Parameters
    pgg = [pgg0, zeros(1,nPoints-1)];
    pgi = [pgi0, zeros(1,nPoints-1)];
    mgg = [mgg0, zeros(1,nPoints-1)];
    mgi = [mgi0, zeros(1,nPoints-1)];
    Tgg = [Tgg0, zeros(1,nPoints-1)];
    Tgi = [Tgi0, zeros(1,nPoints-1)];

    for m = 1:nPoints-1

        % SET OLD VALUES
        dub_old = dub;
        pw0_old = pw0;
        rhow0_old = rhow0;

        % TIME-DIFFERENTIAL QUANTITIES
        % Temperature Change from Joule-Thomson Effect
        if ~ismember('jt',excludePhysics)
            DTjtg = kjtg*(pgg0 - pf0); % temperature reduction of air mass transferred to bubble (generator) [K] NOTE: DTg = dTg*Dt
            DTjti = kjti*(pgi0 - pf0); % temperature reduction of air mass transferred to bubble (injector) [K] NOTE: DTg = dTg*Dt
            Tb0 = (Tb0*mb0 + (Tgg0 - DTjtg)*dmag*Dt + (Tgi0 - DTjti)*dmai*Dt)...
                /(mb0 + dmag*Dt + dmai*Dt); % bubble temperature after differentials of air mass enter it [K]
        end

        % Bubble Depth
        if ~ismember('zb',excludePhysics)
            dzb = 2*g*rbcu*Dt/rb0.^3; % differential variation of bubble depth with time [m s-1] 
            rbcu = rbcu + rb0^3; % faster than dzb = 2*g*sum(rb(1:m).^3)*Dt/rb0.^3; 

            % Bubble Detachment from Gun
            if zb0 < z - rg - rb0
                Vgun = 0; % results in smooth detachment
                % rb0 = (3*(Vb0+Vgun)/(4*pi))^(1/3); % add to Vgun = 0 for abrup detachment
            end
        else
            dzb = 0;
        end

        % Bubble Volume    
        dVb = 4*pi*ub0*rb0^2;  % differential variation of volume with time [m3 s-1]

        % Mass Transfer of Dry Air
        if ~ismember('ma',excludePhysics)
            % Mass of Dry Air in the Bubble (Generator) [Gained]
            if mag0 <= eta*mgtg && pgg0 > pf0 % before all possible air mass in the generator chamber is released (mgtg*(1-eta) remains in the chamber)
                dmag = taug*(Vgg^beta) * sqrt((mgg0/Vgg)*(pgg0-pf0)); % [Kg s-1]
%                 taug = 2e-3; % PROVISIONAL !!!
%                 dmag = mgg0/taug; % PROVISIONAL !!!
            else % after all possible air mass is released (mgt*(1-eta) remains in the chamber)
                dmag = 0;
            end 

            % Mass of Dry Air in the Bubble (Injector) [Gained]
            if m > round(ti*fsp+1) && mai0 <= eta*mgti && pgi0 > pf0  % after injection time ti and before all possible air mass in the injector chamber is released (mgti*(1-eta) remains in the chamber)
                dmai = taui*(Vgi^beta) * sqrt((mgi0/Vgi)*(pgi0-pf0)); % [Kg s-1]
%                 taui = 40e-3;  % PROVISIONAL !!!
%                 dmai = mgi0/taui;  % PROVISIONAL !!!
            else % after all possible air mass is released (mgt*(1-eta) remains in the chamber)
                dmai = 0;
            end 
        else
            dmag = 0;
            dmai = 0;
        end

        % Mass of Water Vapour in the Bubble [Gained]
        if ~ismember('mv',excludePhysics)
            dmv = kec*krou*4*pi*rb0^2/sqrt(2*pi*Rv)*(pveq/sqrt(Tw) ...
                - kcor*pv0/sqrt(Tb0)); % [Kg s-1]
        else
            dmv = 0;
        end

        % Total Mass in the Bubble
        dma = dmag + dmai;
        dmb = dma + dmv; % [Kg s-1]

        % Heat transferred into Water through Bubble Wall [Lost] 
        if ~ismember('qw',excludePhysics)
            dQw = epsil*4*pi*rb0^2*(Tb0-Tw); % differential variation of heat with time [J s-1]
        else
            dQw = 0;        
        end

        % Heat Contained in the Evaporated Water [Lost]
        if ~ismember('qv',excludePhysics)
            dmvin = dmv*(sign(dmv)+1)/2; % dmvin = dmv for dmv > 0, and dmvin = 0 for dmv < 0 (i.e. heat is not recovered when vapour condenses)
            dQv = dmvin*Lv; % [J s-1]
        else
            dQv = 0;
        end

        % Total Heat [Lost]
        dQb = dQv + dQw; % [J s-1]

        % Temperature Inside the Bubble
        Rb = (ma0*Ra + mv0*Rv)/(ma0 + mv0); % specific gas constant of gas mixture [J kg-1 K-1]
        dTb = (Rb*Tb0*dmb -dQb -pb0*dVb)/(ma0*cva + mv0*cvv); % differential variation of temperature with time [K s-1]

        % Specific Enthalpy of Water at the Bubble Wall
        % [similar to dhb = (hb0 - hb0_old)/Dt but accurate only if Dt large]
        dpa = Ra*Tb0/Vb0*dma + Ra*ma0/Vb0*dTb - Ra*ma0*Tb0/Vb0^2*dVb;
        dpv = Rv*Tb0/Vb0*dmv + Rv*mv0/Vb0*dTb - Rv*mv0*Tb0/Vb0^2*dVb;
        dpf = dpa + dpv; % NOTE: alternative dpf = Rb*Tb0/Vb0*dmb + Rb*mb0/Vb0*dTb - Rb*mb0*Tb0/Vb0^2*dVb
        dpb = dpf + ub0/rb0^2*(2*sigmaw + 4*muwe*ub0) - 4*muwe/rb0*dub; % pressure outside bubble wall [Pa s-1]
        dhb = dpb/rhow0 - hb0*drhow/rhow0 - dpw/rhow0; % specific enthalpy [J kg-1 s-1]

        % Radial Velocity and Acceleration of the Bubble
        dub = bubbleAcceleration(rb0,ub0,pb0,hb0,dpb,dhb,cw0,rhow0,pw0); % particle acceleration in the bubble wall [m s-2]
        d2ub = (dub - dub_old)/Dt; % time derivative of the particle acceleration in the bubble wall [m s-3]

        % NEXT VALUES     
        % Bubble Parameters
        ub0 = ub0 + dub*Dt + d2ub*Dt^2/2; % next bubble wall velocity [m s-1] NOTE!: 2nd derivative slightly reduces the gap between peaks and valleys
        rb0 = rb0 + ub0*Dt + dub*Dt^2/2; % next bubble radius [m] NOTE!: 2nd derivative slightly reduces the gap between peaks and valleys
        Tb0 = Tb0 + dTb*Dt; % next temperature in the bubble [K]
        mb0 = mb0 + dmb*Dt; % next mass of gas mixture in the bubble [kg]
        zb0 = zb0 - dzb*Dt; % next bubble depth [m]
        Vb0 = 4/3*pi*rb0^3 - Vgun; % next bubble volume [m3]  

        % Bubble Mass
        mag0 = mag0 + dmag*Dt; % next mass of dry air in the bubble (from generator) [kg]
        mai0 = mai0 + dmai*Dt; % next mass of dry air in the bubble (from injector) [kg]
        ma0 = ma0 + dma*Dt; % next mass of dry air in the bubble (ma0 = mag0 + mai0) [kg]
        mv0 = mv0 + dmv*Dt; % next mass of water vapour in the bubble [kg]

        % Bubble Pressure
        pag0 = mag0*Tb0*Ra/Vb0; % partial pressure of non-condensable gas (dry air) inside the bubble (from generator) [Pa]
        pai0 = mai0*Tb0*Ra/Vb0; % partial pressure of non-condensable gas (dry air) inside the bubble (from injector) [Pa]
        pa0 = ma0*Tb0*Ra/Vb0; % partial pressure of non-condensable gas (dry air) inside the bubble (pa0 = pag0 + pai0) [Pa]
        pv0 = mv0*Tb0*Rv/Vb0; % partial pressure of condensable gas (water vapour) inside the bubble [Pa]
        pf0 = pa0 + pv0; % pressure of gas mixture (dry air + water vapour) insite the bubble [Pa]
        pten = 2*sigmaw/rb0; % pressure due to tension at bubble wall [Pa]
        pvis = 4*muwe*ub0/rb0; % pressure due to viscosity of water at bubble wall [Pa]
        pb0 = pf0 - pten - pvis; % pressure at the external face of the bubble wall [Pa]

        % Water Parameters
        pw0 = patm + rhow_inf*g*zb0; % next hydrostatic pressure [Pa] NOTE: using pw0 (time dep.) introduces a small error, since dhb is calculated assumming pw0 = cst
        pe0 = pw0 + ps(m+1); % next effective ambient pressure [Pa]
        rhow0 = rhow_inf*(((pb0+B)/(pe0+B))^(1/n)); % perturbed water density [Kg m-3]
        cw0 = cw_inf*((pb0+B)/(pe0+B))^((n-1)/(2*n)); % perturbed water sound speed [m s-1]
        hb0 = (pe0+B)/rhow_inf * n/(n-1) * (((pb0+B)/(pe0+B))^((n-1)/n) -1);  % next specific enthalpy (= enthalpy / mass) [J kg-1]

        % Effective Viscosity
        if ~ismember('mu',excludePhysics)
            Re = rhow0*abs(ub0)*2*rb0/muw; % Reynolds number
            muwe = muw*(1 + kvis*Re); % effective dynamic viscosity of water [W m-1 K-1]
        else
            muwe = muw;
        end

        % Gun Parameters
        % [NOTE: Tgg0 and Tgi0 assumed constant]   
        mgg0 = mgtg - mag0; % next mass in the gun chamber (generator) [kg]
        mgi0 = mgti - mai0; % next mass in the gun chamber (injector) [kg]
        pgg0 = mgg0*Ra*Tgg0/Vgg;
        pgi0 = mgi0*Ra*Tgi0/Vgi;

        % Other Parameters
        kjtg = Vgg/(mgg0*cpa)*(kte*Tgg0 - 1); % Joule-Thomson coefficient (generator) [K Pa-1]
        kjti = Vgi/(mgi0*cpa)*(kte*Tgi0 - 1); % Joule-Thomson coefficient (injector) [K Pa-1]
        dpw = (pw0 - pw0_old)/Dt;
        drhow = (rhow0 - rhow0_old)/Dt;

        % Correct Parameters if Gun Volumes Set to Zero
        if Vgg == 0, pgg0 = 0; mgg0 = 0; kjtg = 0; end
        if Vgi == 0, pgi0 = 0; mgi0 = 0; kjti = 0; end

        % BUILD OUTPUT VECTOR (next value)
        % Water Parameters
        rhow(m+1) = rhow0;
        cw(m+1) = cw0;
        pw(m+1) = pw0; 
        pe(m+1) = pe0;

        % Bubble Parameters
        pb(m+1) = pb0;
        ub(m+1) = ub0;
        rb(m+1) = rb0;
        Tb(m+1) = Tb0;
        mb(m+1) = mb0;
        hb(m+1) = hb0;
        Vb(m+1) = Vb0;
        zb(m+1) = zb0;

        % Bubble Pressure
        pf(m+1) = pf0;
        pa(m+1) = pa0;
        pv(m+1) = pv0;
        pag(m+1) = pag0;
        pai(m+1) = pai0;

        % Bubble Mass
        ma(m+1) = ma0;
        mv(m+1) = mv0;
        mag(m+1) = mag0;
        mai(m+1) = mai0;

        % Gun Parameters
        pgg(m+1) = pgg0;
        pgi(m+1) = pgi0;
        mgg(m+1) = mgg0;
        mgi(m+1) = mgi0;
        Tgg(m+1) = Tgg0;
        Tgi(m+1) = Tgi0;
    end

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

    % POPULATE 'Bubble' STRUCTURE
    Bubble.Config = struct(...
        'gunIds',Config.gunIds,...
        'duration',Config.duration,...
        'sampleRate',Config.sampleRate,...
        'waterDensity',Config.waterDensity,...
        'waterSoundSpeed',Config.waterSoundSpeed,...
        'waterTemperature',Config.waterTemperature,...
        'waterSalinity',Config.waterSalinity,...
        'bubbleMethod',Config.bubbleMethod,...
        'bubbleFormula',Config.bubbleFormula,...
        'interactionGhost',Config.interactionGhost);

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
