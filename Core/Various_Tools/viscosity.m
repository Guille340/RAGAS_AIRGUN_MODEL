
%  mu = viscosity(gas,T,varargin)
%
%  DESCRIPTION: Calculates the dynamic viscosity of air or water vapour
%
%  INPUT VARIABLES
% - gas: string specifying the gas type. Two options:
%   ¬ 'air': dry air
%   ¬ 'h2o': water vapour
% - T: temperature of the gas [K]
% - rhov (varargin{1}): density of water vapour [kg m-3]
%
%  OUTPUT VARIABLES
%  - mu: dynamic viscosity of air (gas = 'air') or water vapour (gas = 
%    'h2o') [Pa s]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  FUNCTION CALLS
%  1) mu = viscosity('air',T)
%  2) mu = viscosity('h2o',T,rhov) 
%
%  CONSIDERATIONS & LIMITATIONS
%  - The viscosity of air is temperature dependent.
%  - The viscosity of steam is temperature and pressure dependent.
%  - The dependence on pressure is controlled by the density rhov.
%
%  REFERENCES
%  - Steam properties calculator -> http://www.peacesoftware.de/
%    einigewerte/calc_dampf.php5
%  - Sato and Minamiyama, "Viscosity of Steam at High Temperatures and
%    Pressures" (1964)

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  26 Sep 2016

function mu = viscosity(gas,T,varargin)

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
    case 'air' % Southerland equation
        C = 120; % Sutherland constant [K]
        T0 = 291.15; % reference temperature [K]
        mu0 = 1.827*1e-5; % reference viscosity [Pa s]
        mu = mu0*(T0 + C)/(T + C)*(T/T0)^1.5;
    case 'h2o' % from Sato & Minamiyama (1964)
        g = 9.80665; % gravity force [m s-2]
        mu_atm = 1.751*T^1.11*g*1e-9;
        a = 5.18*10^-7;
        b = -2.78*10^-9;
        c = 9.44*10^-9;
        gam = rhov*g;
        mu = mu_atm + a*(gam/100) + b*(gam/100)^2 + c*(gam/100)^3;
end

        