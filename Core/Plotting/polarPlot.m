% h = polarPlot(angles,radii,varargin)
%
% DESCRIPTION
% Polar plot of the curve defined by the input ANGLES and RADII. The radial
% limits, the legend strings, or the appearance of the polar plot itself
% can be modified through available variable input arguments.
%
% INPUT ARGUMENTS
% - angles: vector of angles in mathematical format ("0 deg E, 90 deg N") [deg]
%   If ANGLES is a matrix, each column represents a curve.
% - radii: vector of radii of same size as ANGLES. If RADII is a matrix, each
%   column represents a curve.
% - radiusLim: two-element vector specifying the radial limits.
% - legendText: cell vector contaning the name of the polar curves.
% - nCircles: number of concentric circles in the polar plot (default = 5)
% - nRadials: number of radial transects in the polar plot (default = 13)
%
% OUTPUT ARGUMENTS
% - h: handle to polar plot
%
% FUNCTION CALL
% 1. h = polarPlot(angles,radii)
% 2. h = polarPlot(angles,radii,radiusLim)
% 3. h = polarPlot(angles,radii,radiusLim,legendText)
% 4. h = polarPlot(angles,radii,radiusLim,legendText,nCircles)
% 5. h = polarPlot(angles,radii,radiusLim,legendText,nCircles,nRadials)
% 
% FUNCTION DEPENDENCIES
% - None
%
% TOOLBOX DEPENDENCIES
% - MATLAB (Core)
%
% See also SEMIPOLARPLOT, PLOTHORIZONTALDIRECTIVITY

% VERSION 1.0
% Author: Guillermo Jimenez-Arranz
% Date: 08 Jun 2022
% Email: gjarranz@gmail.com

function h = polarPlot(angles,radii,varargin)

% ANGLES: mathematical format (0 E, 90 N)

narginchk(0,6)

switch nargin
    case 2
        radiusLim = [0 1];
        legendText = [];
        nCircles = 5;
        nRadials = 13;
    case 3
        radiusLim = varargin{1};
        legendText = [];
        nCircles = 5;
        nRadials = 13;
    case 4
        radiusLim = varargin{1};
        legendText = varargin{2};
        nCircles = 5;
        nRadials = 13;
    case 5
        radiusLim = varargin{1};
        legendText = varargin{2};
        nCircles = varargin{3};
        nRadials = 13;   
    case 6
        radii = varargin{1};
        radiusLim = varargin{2};
        legendText = varargin{3};
        nCircles = varargin{4};
        nRadials = 13;
end

% Make Radius a Column Vector
if isvector(angles)
    angles = angles(:);    
end
if isvector(radii)
    radii = radii(:);    
end

% General Parameters
rmin = min(radiusLim);
rmax = max(radiusLim);
Dr = rmax - rmin;
ang_direc = angles*pi/180; % directivity angles [rad]

% Initialise Figure
h = figure;
hold on

% Plot Circle Polygon
nPoints = 72;
ang_circles = (3*pi/2:-2*pi/(nPoints-1):-pi/2)';
x_pol = Dr.*cos(ang_circles);
y_pol = Dr.*sin(ang_circles);
fill(x_pol,y_pol,[1 1 1],'LineWidth',1.5)

% Plot Radials
radialStep = 2*pi/(nRadials-1);
ang_radials = (3*pi/2:-radialStep:-pi/2)';
theta_radials = (-pi:radialStep:pi)';
for m = 1:nRadials
    x_radial = Dr*cos(ang_radials(m));
    y_radial = Dr*sin(ang_radials(m));
    x_text = 1.08*x_radial - 0.03*Dr;
    y_text = 1.08*y_radial + 0.00*Dr;
    plot([0 x_radial],[0 y_radial],'Color',[0.7 0.7 0.7])
    htxt = text(x_text,y_text,sprintf('%0.0f°',theta_radials(m)*180/pi));
end
htxt.String = ''; % remove last text-type angle

% Plot Concentric Circles
r_step = (rmax - rmin)/(nCircles - 1);
r_circles = rmin:r_step:rmax;
for m = 1:nCircles
    x_circle = r_circles(m)*cos(ang_circles);
    y_circle = r_circles(m)*sin(ang_circles);
    x_text = 0.02;
    y_text = -r_circles(m) - 0.02;
    plot(x_circle,y_circle,'Color',[0.7 0.7 0.7])
    text(x_text,y_text,sprintf('%0.2f',r_circles(m)),'Color',[0.7 0.7 0.7])
end

% External Lines
plot(x_pol,y_pol,'k','LineWidth',1.5)

% Plot Directivity
x_direc = radii.*cos(ang_direc);
y_direc = radii.*sin(ang_direc);
hplo = plot(x_direc,y_direc,'LineWidth',1.5);

% Plot Legend
if ~isempty(legendText)
    hleg = legend(hplo,legendText,'Location','SouthWest');
end
hleg.Position(1) = hleg.Position(1) - 5e-3;

% Adjust
set(h,'color',[1 1 1])
set(gca,'xcolor',[1 1 1],'ycolor',[1 1 1])
axis([-1 1 -1 1])
pbaspect([1 1 1])
set(gcf,'units','normalized','outerposition',[0.1 0.05 0.8 0.9])