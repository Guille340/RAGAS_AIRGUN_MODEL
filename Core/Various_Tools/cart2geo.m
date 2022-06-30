% [dist,sigma,theta] = cart2geo(x,y,z)
%
% DESCRIPTION
% Calculates the distance DIST, elevation THETA, and azimuth SIGMA
% relative to the origin of coordinated (0,0,0) from a point determined 
% by its the cartesian coordinates (X,Y,Z).
%
% Use this function to calculate the following parameters from the
% PLOT.JSON file: "souToRecDistance" (DIST), "souToRecAzimuth" (SIGMA)
% and "souToRecElevation" (THETA).

function [dist,sigma,theta] = cart2geo(x,y,z)

dist = sqrt(x.^2 + y.^2 + z.^2);
sigma = atan2(x,y) * 180/pi;
theta = asin(z./dist) * 180/pi;

