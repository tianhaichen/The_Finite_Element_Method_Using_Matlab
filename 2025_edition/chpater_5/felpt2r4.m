function [m] = felpt2r4(xleng,yleng)

%-------------------------------------------------------------------
%  Purpose:
%     element matrix for transient term of two-dimensional 
%     Laplace's equation using linear triangular element
%
%  Synopsis:
%     [m]=felpt2r4(xleng,yleng)
%
%  Variable Description:
%     m - element stiffness matrix (size of 3x3)   
%     xleng - element size in the x-axis
%	  yleng - element size 	in the y-axis 
%-------------------------------------------------------------------

% element matrix

 A=xleng*yleng; % area of the triangule
 
 m = (A/36)* [4 2 1 2;
              2 4 2 1;
              2 2 4 2;
			  2 1 2 4];