% muw = dynamicViscosityWater(T,S)
%
% DESCRIPTION
% Calculates the absolute (dynamic) viscosity of sea water MUW at 1 atm 
% based on an input temperature T and salinity S.
%
% INPUT ARGUMENTS
% T: sea water temperature [K]
% S: sea water salinity [ppt]
%
% OUTPUT ARGUMENTS
% muw: dynamic viscosity of sea water [kg m-1 s-1]
%
% FUNCTION CALL
% muw = dynamicViscosityWater(T,S)
% 
% FUNCTION DEPENDENCIES
% - None
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% REFERENCES
% Sharqawy, M.H., Lienhard, J.H., and Zubair, S.M. (2010) "Thermophysical 
% properties of seawater: A review of existing correlations and data", 
% Desalination and Water Treatment, 16(1-3), 354-380.
% https://doi.org/10.5004/dwt.2010.1079
%
% See also SURFACETENSIONWATER

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function muw = dynamicViscosityWater(T,S)

% Unit Conversion
T = T - 273.15; % sea water temperature [C]
S = S * 1e-3; % sea water salinity [kg/kg]

% Dynamic Viscosity of Pure Water at 1 atm [kg m-1 s-1]
mup = 4.2844e-5 + (0.157*(T + 64.993).^2 - 91.296).^-1;

% Surface Tension of Sea Water at 1 atm [kg m-1 s-1]
A = 1.541 + 1.998e-2*T - 9.52e-5*T.^2;
B = 7.974 - 7.561e-2*T + 4.724e-4*T.^2;
muw = mup.*(1 + A.*S + B.*S.^2);
