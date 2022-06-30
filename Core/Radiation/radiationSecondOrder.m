% [p,u,t] = radiationSecondOrder(Dynamics,distance,sampleRate)
%
% DESCRIPTION
% Computes the acoustic signature of the individual air gun at the specified 
% DISTANCE. The motion parameters of the air gun bubble are given in the 
% DYNAMICS structure. The function uses the second order approximation from
% Gilmore (1952). In this formula, the SAMPLERATE is used to improve avoid
% nulls by adding second differential terms.
%
% INPUT ARGUMENTS
% - Dynamics: single-element structure containing the bubble motion parameters
%   for a single air gun. The structure is obtained from the RADIATION
%   structure produced by BUBBLERADIATIONS.m. For details about the content
%   of this structure see INITIALISEBUBBLESTRUCT.
% - distance: distance, in metres, between the centre of the array and the
%   receiver.
% - sampleRate: sampling rate for the time-dependent bubble motion parameters,
%   in Hz.
%
% OUTPUT ARGUMENTS
% - p: time-dependent pressure signature of the air gun bubble [Pa]
% - u: time-dependent velocity signature of the air gun bubble [m s-1]
% - t: time axis for the acoustic signatures [s]
%
% FUNCTION CALL
% 1. [p,u,t] = radiationSecondOrder(Dynamics,distance)
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
% See also BUBBLERADIATION, RADIATIONFIRSTORDER

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function [p,u,t] = radiationSecondOrder(Dynamics,distance,sampleRate)

% Import Bubble Parameters
hb = Dynamics.bubbleEnthalpy;
ub = Dynamics.bubbleVelocity;
rb = Dynamics.bubbleRadius;
cw = Dynamics.waterSoundSpeed;
rhow = Dynamics.waterDensity;
Dt = 1/sampleRate;

% Pressure and Particle Velocity at Receiver
y = rb.*(hb + ub.^2/2);
d2y0 = diff(diff(y));
d2y = [d2y0(1) d2y0(1) d2y0];
y = y + d2y*Dt^2/2; % 2nd derivative added to avoid zero values in y (-> NaN in p and u)
K3 = cw.^3.*rb.^2.*ub./y.^2.*(1-ub.^2./(2*cw.^2)) - cw.^2.*rb./y.*(1 - ub./cw);
u = y./(cw*distance) + K3.*y.^2./(cw.^3*distance^2).*(1 - y./(cw.^2*distance) ...
    + K3.^2.*y.^4./(2*cw.^8*distance^4));
p = rhow.*(y/distance - u.^2/2) + rhow./(2*cw.^2).*(y/distance - u.^2/2).^2;
t = distance./mean(cw,'omitnan');
