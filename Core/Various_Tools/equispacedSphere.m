function [x,y,z] = equispacedSphere(nAzim,radius,origin,nAzimAligned,semisphere)

% N: number of points in the main (largest) circumference

narginchk(1,5)

% Initialise Variables 
switch nargin
    case 1
        radius = 1;
        origin = [0 0 0];
        nAzimAligned = 1;
        semisphere = false;
    case 2
        origin = [0 0 0];
        nAzimAligned = 1;
        semisphere = false;
    case 3
        nAzimAligned = 1;
        semisphere = false;
    case 4
        semisphere = false;
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

% Error Control ('nAzimuthsAligned)
if ~isnumeric(nAzimAligned) || ~isscalar(nAzimAligned) ...
        || nAzimAligned < 1 || nAzimAligned > nAzim
    nAzimAligned = 1;
    warning(['NAZIMUTHSALGINED must be a scalar betwee 1 and NAZIMUTHS. '...
        'NAZIMUTHSALIGNED = %d will be used'],nAzimAligned)
end

% Process Cartesian Coordinates
nElev = floor(nAzim/4); % number of elevation angles
theta_step = pi/(2*nElev); % elevation angle step [rad]
cnt = 1;
for m = 1:nElev
    theta0 = theta_step*(m - 1); % elevation angle [rad]
    z0 = radius*sin(theta0); % depth [m]
    nAzim0 = ceil(nAzim*cos(theta0)/nAzimAligned)*nAzimAligned; % number of azimuthal points in circle
    for n = 1:nAzim0
        azim0 = 2*pi/nAzim0*(n - 1); % azimuth [deg]
        x0 = radius*cos(theta0)*cos(azim0); % x position [m]
        y0 = radius*cos(theta0)*sin(azim0); % y positioni [m]
        
        % Populate Coordinate Vectors
        x(cnt) = x0;
        y(cnt) = y0;
        z(cnt) = z0;
        cnt = cnt + 1;
    end
end

% Last Point
theta0 = theta_step*nElev; % elevation angle [rad]
sigma0 = 0; % azimuth [deg]
x(cnt) = radius*cos(theta0)*cos(sigma0); % x position [m]
y(cnt) = radius*cos(theta0)*sin(sigma0); % y positioni [m]
z(cnt) = radius*sin(theta0); % depth [m]

% Convert into Columns
x = x(:);
y = y(:);
z = z(:);

% Complete Sphere
if ~semisphere
    x_mirror = x(nAzim + 1:end);
    y_mirror = y(nAzim + 1:end);
    z_mirror = -z(nAzim + 1:end);
    
    x = [x; x_mirror];
    y = [y; y_mirror];
    z = [z; z_mirror];
end

% Shift to Origin
x = x + origin(1);
y = y + origin(2);
z = z + origin(3);
