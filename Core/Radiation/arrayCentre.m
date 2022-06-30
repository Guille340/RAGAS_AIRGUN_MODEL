% [xc,yc,zc] = arrayCentre(Bubble)
%
% DESCRIPTION
% Returns the central (mean) position of the air gun array defined in the
% BUBBLE structure.
%
% INPUT ARGUMENTS
% - Bubble: structure containing the input configuration information and 
%   predicted time-depended bubble motion parameters. The structure is 
%   generated with functions MULTIPHYSICSDYNAMICS.m or POLYTROPICDYNAMICS.m.
%
% OUTPUT ARGUMENTS
% - xc: x coordinate of the centre of the array [m]
% - yc: y coordinate of the centre of the array [m]
% - zc: z coordinate of the centre of the array [m]
%
% FUNCTION CALL
% 1. [xc,yc,zc] = arrayCentre(Bubble)
% 
% FUNCTION DEPENDENCIES
% - None
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% See also BUBBLERADIATION

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function [xc,yc,zc] = arrayCentre(Bubble)

xc = mean([Bubble.Array.gunX]);
yc = mean([Bubble.Array.gunY]);
zc = mean([Bubble.Array.gunZ]);
