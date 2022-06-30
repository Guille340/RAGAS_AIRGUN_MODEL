%  y = COMPACTNUM(x,varargin)
%
%  DESCRIPTION
%  Converts the input numerical vector X into a matrix of alpha-numerical 
%  strings. Each string consists of the  number expressed as floating point 
%  with NDEC decimal places followed by a metric prefix letter. For example, 
%  for X = 2525 and NDEC = 2, Y = '2.53k'. This function is particularly useful 
%  for showing in the x-axis the central frequencies of a band spectrum in a 
%  neater way.
% 
%  INPUT VARIABLES
%  - x: vector of decimal numbers.
%  - ndec (varargin{1}): number of decimal places to show in the strings.
%   
%  OUTPUT VARIABLES
%  - y: alpha-numerical matrix of strings showing the compact form of input
%    numbers X.
%
%  FUNCTION CALL
%  y = COMPACTNUM(X)
%  y = COMPACTNUM(X,NDEC)
%
%  NDEC can be left empty ([]) or omitted from the call, its default 
%  value will be used (NDEC = 2).
%
%  FUNCTION DEPENDENCIES
%  - None
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  EXAMPLE
%  x = [0.001 0.01 0.1 1 11 111 1111 11111 111111 1111111]
%  y = compactnum(x,2)

%  VERSION 1.1
%  Revised by: Guillermo Jimenez Arranz
%  Date: 20 May 2021
%  - Added help
%  - Improved error control
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  30 Jan 2018

function y = compactnum(x,varargin)

% Error Control
narginchk(0,3)
ndec = 2; % default value
if nargin == 2
    ndec = varargin{1}; % maximum number of decimal digits
end

if ~isvector(x) || ~isnumeric(x)
    error('Input X must be a numerical vector')
end

if isempty(ndec)
    ndec = 2;
end

% General Parameters
x = x(:); % convert input into column vector
npow = floor(floor(log10(abs(x)))/3)*3;
npow(~x) = 0;

% Numerical Prefixes of International System (SI)
prefix = strings(length(x),1);
prefix(npow == -24) = 'y';
prefix(npow == -21) = 'z';
prefix(npow == -18) = 'a';
prefix(npow == -15) = 'f';
prefix(npow == -12) = 'p';
prefix(npow == -9) = 'n';
prefix(npow == -6) = '\mu';
prefix(npow == -3) = 'm';
prefix(npow == 0) = '';
prefix(npow == 3) = 'k';
prefix(npow == 6) = 'M';
prefix(npow == 9) = 'G';
prefix(npow == 12) = 'T';
prefix(npow == 15) = 'P';
prefix(npow == 18) = 'E';
prefix(npow == 21) = 'Z';
prefix(npow == 24) = 'Y';
prefix = cellstr(prefix);

compnum = round(x./10.^npow * 10^ndec)/10^ndec;
compstr = strtrim(cellstr(num2str(compnum)));
y = strcat(compstr,prefix);
