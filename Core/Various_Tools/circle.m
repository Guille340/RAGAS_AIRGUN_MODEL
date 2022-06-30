
function [x,y] = circle(x0,y0,r,varargin)

%*************************************************************************
%  [x,y] = circle(x0,y0,r,varargin)
%  
%  DESCRIPTION: returns the horizontal and vertical position (¦x¦,¦y¦)
%  for each point on a circle with origin at (¦x0¦,¦y0¦) and radius ¦r¦.
%  If ¦r¦ is a vector, each column in ¦x¦ and ¦y¦ corresponds to a 
%  different circle. The number of points is by default 361 (1 degree 
%  resolution)
%
%  INPUT VARIABLES
%  - x0: horizontal position of the origin of coordinates [m]
%  - y0: vertical position of the origin of coordinates [m]
%  - r: radius or vector of radii [m]
%  - M (varargin{1}): number of points in the circle. By default, M = 361
%
%  OUTPUT VARIABLES
%  - x: vector of horizontal positions in the circle [m]. If ¦r¦ is a
%    vector, ¦x¦ is an array with as many columns as elements in ¦r¦ (one
%    column per circle)
%  - y: vector of vertical positions in the circle [m]. If ¦r¦ is a
%    vector, ¦y¦ is an array with as many columns as elements in ¦r¦ (one 
%    column per circle)
%
%  FUNCTION CALLS
%  1) [x,y] = circle(x0,y0,r)
%     ¬ M = 361
%  2) [x,y] = circle(x0,y0,r,M)
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  14 Jun 2016
%
%************************************************************************* 

switch nargin
    case {0 1 2}
        error('Not enough input arguments')
    case 3
        M = 361; % number of points (L = 360/ares + 1, ares = angle resolution)
    case 4
        M = varargin{1};
    otherwise
        error('Too many output arguments')
end

r = r(:)'; % row vector
N = length(r);
ang = (0:360/(M-1):360)'*pi/180;
r = repmat(r,M,1);
ang = repmat(ang,1,N);

x = x0 + r.*cos(ang);
y = y0 + r.*sin(ang);
