% sigmaw = surfaceTensionWater(T,S)
%
% DESCRIPTION
% Calculates the surface tension of sea water SIGMAW at 1 atm based on
% an input temperature T and salinity S.
%
% INPUT ARGUMENTS
% T: sea water temperature [K]
% S: sea water salinity [ppt]
%
% OUTPUT ARGUMENTS
% sigmaw: surface tension of sea water [N m-1]
%
% FUNCTION CALL
% sigmaw = surfaceTensionWater(T,S)
%
% FUNCTION DEPENDENCIES
% - None
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% REFERENCES
% Nayar, K. G., Panchanathan, D., McKinley, G. H., and Lienhard, J.H.
% "Surface tension of seawater", Journal of Physical and Chemical Reference
% Data, 43(4), https://doi.org/10.1063/1.4899037
%
% See also DYNAMICVISCOSITYWATER

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022

function sigmaw = surfaceTensionWater(T,S)

% Surface Tension of Pure Water at 1 atm [N m-1]
sigmap = 0.2358 * (1 - T/647.096)^1.256 * (1 - 0.625*(1 - T/647.096)); 

% Surface Tension of Sea Water at 1 atm [N m-1]
sigmaw =  sigmap * (1 + 3.766e-4*S + 2.347e-6*S*(T-273.15));
