%----------------------------------------------------------------
% Example 6.5.2
%   Gauss-Legendre quadrature of a function in 2-dimension
%
% Problem description
%   integrate f(x,y) = 1+4xy-3x^2*y^3 + x^4*y^6 over -1<x<1 and -1<y<1
%
% Varaible descriptions 
%   point2 = integration points 
%   weight2 = weighting coefficients 
%   ng1x = number of integration points along x-axis
%   ng1y = number of integration points along y-axis
clc
clear
ng1x = 3;       % (2*ng1x-1) = 4
ng1y = 4;       % (2*ng1y-1) = 6

[point2,weight2]=feglqd2(ng1x,ng1y);        % integration points and weights

%-----------------------------------------
% summation for numerical integration
%-----------------------------------------

value = 0.0;
for intx = 1:ng1x
   x=point2(intx,1);        % sampling point in x-axis
   wtx=weight2(intx,1);     % weight in x-axis
   for inty=1:ng1y
      y=point2(inty,2)
      wty=weight2(inty,2);
      func=1+4*x*y-3*x^2*y^2+x^4*y^6;
      value = value + func*wtx*wty;
   end
end

value   % print the solution
%-------------------------------------------------------------------
