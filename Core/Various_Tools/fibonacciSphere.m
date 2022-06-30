function [x,y,z] = fibonacciSphere(nPoints)

phi = pi*(3 - sqrt(5));

for m = 1:nPoints
    y0 = 1 - 2*(m/(nPoints));
    r = sqrt(1 - y0.^2); % radius
    
    theta = phi*m;
    x0 = r*cos(theta);
    z0 = r*sin(theta);
    
    x(m) = x0;
    y(m) = y0;
    z(m) = z0;  
end