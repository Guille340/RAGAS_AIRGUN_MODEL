% [x,y,z] = geo2cart(dist,sigma,theta)
%
% DESCRIPTION
% Calculates the cartesian coordinates (X,Y,Z) for a point determined
% by its distance DIST, elevation THETA, and azimuth SIGMA from the
% origin.

function [x,y,z] = geo2cart(dist,sigma,theta)

theta = theta * pi/180;
sigma = sigma*pi/180;

x = dist*cos(theta)*sin(sigma);
y = dist*cos(theta)*cos(sigma);
z = dist*sin(theta);