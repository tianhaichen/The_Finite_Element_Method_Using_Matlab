function [k] = felp2dt3(x1,y1,x2,y2,x3,y3)
%---------------------------------------------------------------------
% Purpose:
%	element matrix for two-dimensional Laplace's equation
% 	using three-node linear trianglular element 
% Synopsis:
%	[k]=felp2dt3(x1,y1,x2,y2,x3,y3)
%
% Variable Description:
%     k - element stiffness matrix (size of 3x3)   
%     x1, y1 - x and y coordinate values of the first node of element
%     x2, y2 - x and y coordinate values of the second node of element
%     x3, y3 - x and y coordinate values of the third node of element
%---------------------------------------------------------------------

% element matrix 

A=0.5*(x2*y3+x1*y2+x3*y1-x2*y1-x1*y3-x3*y2); % area of the triangle

k11=((x3-x2)^2+(y2-y3)^2)/(4*A);
k12=((x3-x2)*(x1-x3)+(y2-y3)*(y3-y1))/(4*A);
k13=((x3-x2)*(x2-x1)+(y2-y3)*(y1-y2))/(4*A);
k21=k12;
k22=((x1-x3)^2+(y3-y1)^2)/(4*A);
k23=((x1-x3)*(x2-x1)+(y3-y1)*(y1-y2))/(4*A);
k31=k13;
k32=k23;
k33=((x2-x1)^2+(y1-y2)^2)/(4*A);
k=[k11,k12,k13;
    k21,k22,k23;
    k31,k32,k33];