% [x,y,z] = gridSphere(radius,origin,sigma,theta)
%
% DESCRIPTION
% Calculates the points around a sphere of specified RADIUS and ORIGIN
% for the given azimuthal and elevantion angles (SIGMA, THETA). The angles
% follow the mathematical (as opposed to geographic) convention, that is,
% "0 deg East, +90 deg North" (as opposed to "0 deg North, +90 deg East"). 
%
% INPUT ARGUMENTS
% - radius: radius of the sphere.
% - origin: three-element vector defining the [x y z] location for the
%   central point of the sphere.
% - sigma: vector of azimuths [deg]
% - theta: vector of elevations [deg]
%
% OUTPUT ARGUMENTS
% - x: colum vector representing the x coordinate of the points in the sphere.
% - y: colum vector representing the y coordinate of the points in the sphere.
% - z: colum vector representing the x coordinate of the points in the sphere.
%
% FUNCTION CALL
% 1. [x,y,z] = gridSphere(radius,origin,sigma,theta)
% 
% FUNCTION DEPENDENCIES
% - None
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function [x,y,z] = gridSphere(radius,origin,sigma,theta)

% SIGMA: azimuth vector (0 deg N, 90 deg E)
% THETA: elevation vector (0 E-W, 90 S)

narginchk(1,4)

% Initialise Variables 
switch nargin
    case 0
        radius = 1;
        origin = [0 0 0];
        sigma = 0:45:315;
        theta = 0:10:90;
    case 1
        origin = [0 0 0];
        sigma = 0:45:315;
        theta = 0:10:90;
    case 2
        sigma = 0:45:315;
        theta = 0:10:90;
    case 3
        theta = 0:10:90;
end
        
% Error Control ('radius')
if ~isnumeric(radius) || ~isscalar(radius)
    radius = 1;
    warning(['Non-valid value for input argument RADIUS. RADIUS = %0.0f '...
        'will be used'],radius) 
end

% Error Control ('origin')
if ~isnumeric(origin) || ~isvector(origin) || length(origin) ~= 3
    origin = [0 0 0];
    warning(['Non-valid value for input argument ORIGIN. ORIGIN = '...
        '[%0.0f %0.0f %0.0f] will be used'],origin(1),origin(2),origin(3)) 
end

% Error Control ('sigma')
if ~isnumeric(sigma) || ~isvector(sigma)
    sigma = 0:45:315;
    warning(['SIGMA must be a numeric vector. SIGMA will contain %d '...
        'equally-spaced azimuthal points'],length(sigma))
end

% Error Control ('sigma')
if ~isnumeric(theta) || ~isvector(theta)
    theta = 0:10:90;
    warning(['THETA must be a numeric vector. THETA will contain %d '...
        'equally-spaced elevation points'],length(sigma))
end

% Convert to Radians
sigma = sigma*pi/180;
theta = theta*pi/180;

% Process Cartesian Coordinates
nElev = length(theta);
nAzim = length(sigma);
cnt = 1;
for m = 1:nElev
    theta0 = theta(m); % elevation angle [rad]
    z0 = radius*sin(theta(m)); % depth [m]
    for n = 1:nAzim
        sigma0 = sigma(n); % azimuth [deg]
        x0 = radius*cos(theta0)*sin(sigma0); % x position [m]
        y0 = radius*cos(theta0)*cos(sigma0); % y positioni [m]
        
        % Populate Coordinate Vectors
        x(cnt) = x0;
        y(cnt) = y0;
        z(cnt) = z0;
        cnt = cnt + 1;
    end
end

% Convert into Columns
x = x(:);
y = y(:);
z = z(:);

% Shift to Origin
x = x + origin(1);
y = y + origin(2);
z = z + origin(3);

% Return Unique Points Only
points = unique([round(x*1e3)*1e-3 round(y*1e3)*1e-3 round(z*1e3)*1e-3],...
    'rows','stable');
x = points(:,1);
y = points(:,2);
z = points(:,3);

