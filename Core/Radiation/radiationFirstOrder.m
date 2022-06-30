% [p,u,t] = radiationFirstOrder(Dynamics,distance)
%
% DESCRIPTION
% Computes the acoustic signature of the individual air gun at the specified 
% DISTANCE. The motion parameters of the air gun bubble are given in the 
% DYNAMICS structure. The function uses the first order approximation from
% Gilmore (1952).
%
% INPUT ARGUMENTS
% - Dynamics: single-element structure containing the bubble motion parameters
%   for a single air gun. The structure is obtained from the RADIATION
%   structure produced by BUBBLERADIATIONS.m. For details about the content
%   of this structure see INITIALISEBUBBLESTRUCT.
% - distance: distance, in metres, between the centre of the array and the
%   receiver.
%
% OUTPUT ARGUMENTS
% - p: time-dependent pressure signature of the air gun bubble [Pa]
% - u: time-dependent velocity signature of the air gun bubble [m s-1]
% - t: time axis for the acoustic signatures [s]
%
% FUNCTION CALL
% 1. [p,u,t] = radiationFirstOrder(Dynamics,distance)
% 
% FUNCTION DEPENDENCIES
% - None
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
% See also BUBBLERADIATION, RADIATIONSECONDORDER

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function [p,u,t] = radiationFirstOrder(Dynamics,distance)

% Import Bubble Parameters
hb = Dynamics.bubbleEnthalpy;
ub = Dynamics.bubbleVelocity;
rb = Dynamics.bubbleRadius;
cw = Dynamics.waterSoundSpeed;
rhow = Dynamics.waterDensity;

% Pressure and Particle Velocity at Receiver
f = -rb.^2./cw.*(hb + ub.^2/2 - ub.*cw);
df = rb.*(hb + ub.^2/2);
u = f/distance^2 + df./(distance*cw);
p = rhow*df/distance - 0.5*rhow*u.^2;
t = distance/cw;
