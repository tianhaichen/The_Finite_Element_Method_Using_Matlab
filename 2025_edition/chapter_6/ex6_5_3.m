%----------------------------------------------------------------
% Example 6.5.2
%   Gauss-Legendre quadrature of a function in 3-dimension
%
% Problem description
%   integrate f(x,y) = 1+4*x^2*y^2-3*x^2*z^4 + y^4*z^6 over -1<(x,y,z)<1 
%
% Varaible descriptions 
%   point3 = integration points 
%   weight3 = weighting coefficients 
%   ng1x = number of integration points along x-axis
%   ng1y = number of integration points along y-axis
%   ng1z = number of integration points along z-axis
clc
clear
ng1x = 2;       % (2*ng1x-1) = 2
ng1y = 3;       % (2*ng1y-1) = 4
ng1z = 4;       %  (2*ng1z-1) = 6

[point3,weight3]=feglqd3(ng1x,ng1y,ng1z);        % integration points and weights

%-----------------------------------------
% summation for numerical integration
%-----------------------------------------

value = 0.0;
for intx = 1:ng1x
   x=point3(intx,1);        % sampling point in x-axis
   wtx=weight3(intx,1);     % weight in x-axis
   for inty=1:ng1y
      y = point3(inty,2);
      wty=weight3(inty,2);
      for intz=1:ng1z
         z=point3(intz,3);
         wtz=weight3(intz,3);
         func=1+4*x^2*y^2-3*x^2*z^4 + y^4*z^6;
         value = value + func*wtx*wty*wtz;
      end
   end
end

value   % print the solution
%-------------------------------------------------------------------
