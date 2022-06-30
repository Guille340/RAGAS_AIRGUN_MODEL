function plotCircle(x,y,r,varargin)

switch nargin
    case {0 1 2}
        error('Not enough input arguments')
    case 3
        colorStr='b';
        nang=100;
    case 4
        colorStr=varargin{1};
        nang=100;
    case 5
        colorStr=varargin{1}; 
        nang=varargin{2};
    otherwise
        error('Too many input arguments')
end

ang=0:2*pi/(nang-1):2*pi;
xp=r*cos(ang);
yp=r*sin(ang);

plot(x+xp,y+yp,colorStr)

end