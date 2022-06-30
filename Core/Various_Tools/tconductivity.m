
function k = tconductivity(gas,T,varargin)

%**************************************************************************
%  k = tconductivity(gas,T,varargin)
%
%  DESCRIPTION: Calculates the thermal conductivity of air or water vapour 
%
%  INPUT VARIABLES
%  - gas: string specifying the gas type. Two options:
%    ¬ 'air': dry air
%    ¬ 'h2o': water vapour
%  - T: temperature of the gas [K]
%  - rhov (varargin{1}): density of water vapour [kg m-3]
%
%  OUTPUT VARIABLES
%  - k: thermal conductivity of air (gas = 'air') or water vapour (gas = 
%    'h2o') [W m-1 K-1]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  FUNCTION CALLS
%  1) k = tconductivity('air',T)
%  2) k = tconductivity('h2o',T,rhov) 
%
%  CONSIDERATIONS & LIMITATIONS
%  - The thermal conductivity of air is temperature dependent.
%  - The thermal conductivity of steam is temperature and pressure dependent.
%  - The dependence on pressure is controlled by the density rhov.
%
%  REFERENCES
%  - Steam properties calculator -> http://www.peacesoftware.de/
%    einigewerte/calc_dampf.php5
%  - Keyes and Vines, "The Thermal Conductivity of Steam", 1964
%  - Brain, "Thermal Conductivity of Steam at Atmospheric Pressure", 1969
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  26 Sep 2016
%
%**************************************************************************

switch gas
    case 'air'
        if nargin < 2
            error('Not enough input arguments')
        elseif nargin > 2
            error('Too many input arguments')
        end
    case 'h2o'
        if nargin < 3
            error('Not enough input arguments')
        elseif nargin == 3
            rhov = varargin{1}; % density of water vapour [kg m-3]
        else
            error('Too many input arguments')
        end
    otherwise
        error('Wrong input string for GAS')
end

switch gas
    case 'air'
        k = -3.933*1e-4 + T*1.0184*1e-4 - T^2*4.857*1e-8;
    case 'h2o' % from Keyes & Vines (1964)
        t = T -273.15; % temperature [ºC]
        k_atm = 0.01842*sqrt(T)./(1 + 5485./(T.*10.^(12./T))); % thermal conductivity at atmospheric pressure [W m-1 K-1] (equivalent to [W m-1 ºC-1])
        s = (sign(t-100)+1)/2;
        c = 20.93*1e-5 + 3.673*1e-11*(s.*(t-100)).^3.5; % below T = 100 ºC, c = c(100ºC)
        a = 9.521*1e-3*10.^(2133./(s.*(T-373.15)+373.15)); % below T = 100 ºC, a = a(100ºC)
        k = k_atm + c.*(exp(a*rhov*1e-3) - 1);
end

