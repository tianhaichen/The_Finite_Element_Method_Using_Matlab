%----------------------------------------------------
% to solve the two-dimensional Laplace's equation given
% as u,xx + u,yy = 0 , 0<x<5, 0<y<10
% u(x,0) = 0, u(x,10) = 100*sin(pi*x/10),
% u(0,y) = 0, u,x(5,y) = 0
% using bilinear rectangluar elements 
% 
% Variable description
%	k = element matrix 
% 	f = element vector 
%	kk = system matrix 
%	ff system vector 
%	gcoord = coordinate values of each node 
%	nodes = nodal connectivity of each element 
%	index = a vector containing system dofs associated with each element 
% 	bcdof = a vector containing dofs associated with boundary conditions 
%	bcval = a vector containing boundary condition values associated with the dofs in 'bcdof'
%--------------------------------------------------------------------------

%----------------------------------
% input data for control parameters
%----------------------------------

clear
clc
nel = 16;			% number of elements
nnel = 4;			% number of nodes per element 
ndof =1;			% number of dofs per node 
nnode=25;			% total number of nodes in system 
sdof = nnode*ndof;	% total system dofs 

% l = 5;						% domain length
% h = 10;						% domain height
% l_n = 4;					% l direction element number 
% w_n = 4;					% w direction element number
% l_size = l / l_n;			% l direction element size
% w_size = w / w_n;			% w direction element size
% nnode = (l_n+1)*(w_n+1);	% total number of nodes in system 
%--------------------------------------------------
% input data for nodal coordinate values 
% gcoord(i,j) where i->node no. and j->x or y 
%--------------------------------------------------

gcoord = [0.0,0.0;1.25,0.0;2.5,0.0;3.75,0.0;5.0,0.0;
		0.0,2.5;1.25,2.5;2.5,2.5;3.75,2.5;5.0,2.5;
		0.0,5.0;1.25,5.0;2.5,5.0;3.75,5.0;5.0,5.0;
		0.0,7.5;1.25,7.5;2.5,7.5;3.75,7.5;5.0,7.5;
		0.0,10.0;1.25,10.0;2.5,10.0;3.75,10.0;5.0,10.0];
		
		
%-------------------------------------------------------
% input data for nodal connnectivity for each element 
% nodes(i,j) where i-> element no. and j-> connected nodes 
%--------------------------------------------------------
nodes = [1 2 7 6;2 3 8 7;3 4 9 8;4 5 10 9;
		6 7 12 11;7 8 13 12;8 9 14 13;9 10 15 14;
		11 12 17 16;12 13 18 17;13, 14 19 18;14 15 20 19;
		16 17 22 21;17 18 23 22; 18 19 24 23;19 20 25 24];
		
%------------------------------------------
% input data for boundary conditions 
%-----------------------------------
bcdof=[1 2 3 4 5 6 11 16 21 22 23 24 25];
bcval=[0 0 0 0 0 0 0 0 0 38.2683 70.7107 92.388 100];

%--------------------------------------
% initialization of matrices and vectors 
%--------------------------------------
ff = zeros(sdof,1);			% initialization of system force vector 
kk = zeros(sdof,sdof);		% initialization of system matrix 
index = zeros(nnel*ndof,1);	% initialization of index vector 

%---------------------------------------------------------------
% computation of element matrices and vectors and their assembly
%---------------------------------------------------------------
for iel = 1:nel 		% loop for the total number of elements
	for i=1:nnel 		% loop for number of nodes per element 
		nd(i) = nodes(iel,i);	% 1st connected node for (iel)-th element 
		x(i) = gcoord(nd(i),1); 	% extract x value of the node
		y(i) = gcoord(nd(i),2);		% extract y value of the node
	end
	xleng=x(2)-x(1);
	yleng=y(4)-y(1);
	index = feeldof(nd,nnel,ndof);	% extract system dofs associated with element 
	k = felp2dr4(xleng,yleng);	% compute element matrix 
	kk = feasmbl1(kk,k,index);			% assemble element matrices
end 

%--------------------------
% apply boundary conditions 
%--------------------------
[kk,ff] = feaplyc2(kk,ff,bcdof,bcval);

%----------------------------
% solve the matrix equation 
%----------------------------
fsol = kk\ff;

%----------------------------
% analytical solution 
%----------------------------
for i=1:nnode 
	x=gcoord(i,1);y=gcoord(i,2);
	esol(i) = 100*sinh(0.31415927*y)*sin(0.31415927*x)/sinh(3.1415927);
end 

%-----------------------------
% print both exact and fem solutions
%-----------------------------------

num = 1:1:sdof;
store = [num' fsol esol']