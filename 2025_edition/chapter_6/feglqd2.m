function [point2,weight2]=feglqd2(ng1x,ng1y)
%-----------------------------------------------------
% Purpose:
%   determine the integration points and weighting coefficients
%   of Gauss-Legendre quadrature for two-dimensional integration
%
%   Synopsis:
%       [point2,weight2]=feglqd2(ng1x,ng1y)
%
%   Variable Description
%       ng1x - number of integration points in the x-axis
%       ng1y - number of integration points in the y-axis
%       point2 - vector containing integration points
%       weight2 - vector containing weighting coefficients 
% ------------------------------------------------------------

% determine the largest one between ng1x and ng1y
if ng1x > ng1y
   ng1 = ng1x;
else
    ng1 = ng1y;
end

% initialization

point2 = zeros(ng1,2);
weight2 = zeros(ng1,2);

% find corresponding integration points and weights
[pointx,weightx] = feglqd1(ng1x);
[pointy,weighty] = feglqd1(ng1y);   % quadrature rule for y-axis

% quadrature for two-dimension

for intx = 1:ng1x
   point2(intx,1)=pointx(intx);
   weight2(intx,1)=weightx(intx);
end

for inty = 1:ng1y
   point2(inty,2)=pointy(inty);
   weight2(inty,2)=weighty(inty);
end
