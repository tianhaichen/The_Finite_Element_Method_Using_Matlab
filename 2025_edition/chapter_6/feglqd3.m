function [point3,weight3]=feglqd3(ng1x,ng1y,ng1z)
%-----------------------------------------------------
% Purpose:
%   determine the integration points and weighting coefficients
%   of Gauss-Legendre quadrature for three-dimensional integration
%
%   Synopsis:
%       [point3,weight3]=feglqd3(ng1x,ng1y,ng1z)
%
%   Variable Description
%       ng1x - number of integration points in the x-axis
%       ng1y - number of integration points in the y-axis
%       ng1z - number of integration points in the z-axis
%       point3 - vector containing integration points
%       weight3 - vector containing weighting coefficients 
% ------------------------------------------------------------

% determine the largest one between ng1x and ng1y
if ng1x > ng1y
    if ng1x > ng1z
        ng1 = ng1x;
    else
        ng1=ng1z;
    end
else
    if ng1y > ng1z
        ng1 = ng1y;
    else
        ng1 = ng1z;
    end
end

% initialization

point3 = zeros(ng1,3);
weight3 = zeros(ng1,3);

% find corresponding integration points and weights
[pointx,weightx] = feglqd1(ng1x);
[pointy,weighty] = feglqd1(ng1y);   % quadrature rule for y-axis
[pointz,weightz] = feglqd1(ng1z);   % quadrature rule for z-axis

% quadrature for three-dimension

for intx = 1:ng1x                       % quadrature in x-axis
   point3(intx,1)=pointx(intx);
   weight3(intx,1)=weightx(intx);
end

for inty = 1:ng1y
   point3(inty,2)= pointy(inty);
   weight3(inty,2)= weighty(inty);
end

for intz = 1:ng1z
   point3(intz,3) = pointz(intz);
   weight3(intz,3) = weightz(intz);
end
