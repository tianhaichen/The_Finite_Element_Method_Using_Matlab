%----------------------------------------------------
% to solve the two-dimensional Laplace's equation given
% as u,rr + (u,r)/r + u,zz = 0 , 4<r<6, 0<z<1
% u(4,z) = 100, u,r(6,z) = 20
% u,z(r,0) = 0, u,z(r,1) = 0
% using linear triangular elements 
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
nel = 10;			% number of elements
nnel = 3;			% number of nodes per element 
ndof =1;			% number of dofs per node 
nnode=12;			% total number of nodes in system 
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

gcoord = [4.0,0.0;4.0,1;4.4,0.0;4.4,1.0;
		4.8,0;4.8,1;5.2,0;5.2,1;
		5.6,0;5.6,1;6.0,0;6.0,1];
		
%-------------------------------------------------------
% input data for nodal connnectivity for each element 
% nodes(i,j) where i-> element no. and j-> connected nodes 
%--------------------------------------------------------
nodes = [1 4 2;
		1 3 4;
		3 6 4;
		3 5 6;
		5 8 6;
		5 7 8;
		7 10 8;
		7 9 10;
		9 12 10;
		9 11 12];
		
%------------------------------------------
% input data for boundary conditions 
%-----------------------------------
bcdof=[1 2];
bcval=[100 100];

%--------------------------------------
% initialization of matrices and vectors 
%--------------------------------------
ff = zeros(sdof,1);			% initialization of system force vector 
kk = zeros(sdof,sdof);		% initialization of system matrix 
index = zeros(nnel*ndof,1);	% initialization of index vector 

pi = 4*atan(1);
ff(11)= 120*pi;
ff(12) = 120*pi;
%---------------------------------------------------------------
% computation of element matrices and vectors and their assembly
%---------------------------------------------------------------
for iel = 1:nel 		% loop for the total number of elements
	nd(1) = nodes(iel,1);	% 1st connected node for (iel)-th element 
	nd(2) = nodes(iel,2);	% 2nd connected node for (iel)-th element 
	nd(3) =nodes(iel,3);	% 3rd connected node for (iel)-th element 
	r1 = gcoord(nd(1),1); z1 = gcoord(nd(1),2);	% coord values of 1st node 
	r2 = gcoord(nd(2),1); z2 = gcoord(nd(2),2); % coord values of 2nd node 
	r3 = gcoord(nd(3),1); z3 = gcoord(nd(3),2); % coord values of 3rd node 
	index = feeldof(nd,nnel,ndof);	% extract system dofs associated with element 
	k = felpaxt3(r1,z1,r2,z2,r3,z3);	% compute element matrix 
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
	r=gcoord(i,1);z=gcoord(i,2);
	esol(i) = 100-6*20*log(4) + 6*20*log(r);
end 

%-----------------------------
% print both exact and fem solutions
%-----------------------------------

num = 1:1:sdof;
store = [num' fsol esol']