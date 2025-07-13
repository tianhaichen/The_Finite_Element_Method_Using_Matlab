%----------------------------------------------------
% to solve the transient two-dimensional Laplace's equation given
% as u,xx + u,yy = u,t  , 0<x<5, 0<y<2
% boundary conditions:
% u(0,y,t) = 100, u(5,y,t) = 100,
% u,y(x,0,t) = 0, u,y(x,2,t) = 0
% initial condition:
%  u(x,y,0) = 0 over the domain 
% using bilinear triangular elements and forward difference method
% 
% Variable description
%	k = element matrix for time-independent term(u,xx + u,yy)
%	m = element matrix for time-dependent term (u,t)
% 	f = element vector 
%	kk = system matrix of k 
%	mm = system matrix of m 
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
nel = 8;			% number of elements
nnel = 4;			% number of nodes per element 
ndof =1;			% number of dofs per node 
nnode=15;			% total number of nodes in system 
sdof = nnode*ndof;	% total system dofs 
deltt = 0.1;		% time step size for transient analysis 
stime = 0.0;		% initial time 
ftime = 10;			% termination time 
ntime = fix((ftime-stime)/deltt);

gcoord = [0.0,0.0;1.25,0.0;2.5,0.0;3.75,0.0;5.0,0.0;
		0.0,1.0;1.25,1.0;2.5,1.0;3.75,1.0;5.0,1.0;
		0.0,2.0;1.25,2.0;2.5,2.0;3.75,2.0;5.0,2.0];
		
		
%-------------------------------------------------------
% input data for nodal connectivity for each element 
% nodes(i,j) where i-> element no. and j-> connected nodes 
%--------------------------------------------------------
nodes = [1 2 7 6;2 3 8 7;3 4 9 8;4 5 10 9;
		6 7 12 11;7 8 13 12;8 9 14 13;9 10 15 14];
		
%------------------------------------------
% input data for boundary conditions 
%-----------------------------------
bcdof=[1 5 6 10 11 15];
bcval=[100 100 100 100 100 100];

%--------------------------------------
% initialization of matrices and vectors 
%--------------------------------------
ff = zeros(sdof,1);			% initialization of system force vector 
fn = zeros(sdof,1);			% initialization of effective system vector 
fsol = zeros(sdof,1);		% solution vector 
sol = zeros(2,ntime+1);		% vecotr containing time history solution 
kk = zeros(sdof,sdof);		% initialization of system matrix 
mm = zeros(sdof,sdof);		% initialization of system matrix
index = zeros(nnel*ndof,1);	% initialization of index vector 

%---------------------------------------------------------------
% computation of element matrices and vectors and their assembly
%---------------------------------------------------------------
for iel = 1:nel 		% loop for the total number of elements
	nd(1) = nodes(iel,1);	% 1st connected node for (iel)-th element 
	nd(2) = nodes(iel,2);	% 2nd connected node for (iel)-th element 
	nd(3) = nodes(iel,3);	% 3rd connected node for (iel)-th element 
	nd(4) = nodes(iel,4);	% 4th connected node for (iel)-th element 
	x1 = gcoord(nd(1),1); y1 = gcoord(nd(1),2);	% coord values of 1st node 
	x2 = gcoord(nd(2),1); y2 = gcoord(nd(2),2); % coord values of 2nd node 
	x3 = gcoord(nd(3),1); y3 = gcoord(nd(3),2); % coord values of 3rd node 
	x4 = gcoord(nd(4),1); y4 = gcoord(nd(4),2);	% coord values of 4th node
	xleng = x2 - x1;			% element size in x-axis
	yleng = y4 - y1;			% element size in y-axis
	
	index = feeldof(nd,nnel,ndof);	% extract system dofs associated with element 
	
	k = felp2dr4(xleng,yleng);	% compute element matrix 
	m = felpt2r4(xleng,yleng);     % compute element matrix 
	kk = feasmbl1(kk,k,index);			% assemble element matrices
	mm = feasmbl1(mm,m,index);			% assemble element matrices
end 

%--------------------------
% loop for time integration 
%--------------------------
for in = 1:sdof 
	fsol(in)=0.0;	%initial condition
end

sol(1,1) = fsol(8);						%store time history solution for node no.8
sol(2,1) = fsol(9);						%store time history solution for node no.9
for it = 1:ntime 
	fn = deltt*ff + (mm-deltt*kk)*fsol;	% compute effective column vector 
	[mm,fn] = feaplyc2(mm,fn,bcdof,bcval);	% apply boundary condition
	fsol = mm\fn;
	sol(1,it+1)=fsol(8);
	sol(2,it+1)=fsol(9);
end 

%----------------------------------------
% analytical solution at node 8
%-------------------------------
pi=4*atan(1);
esol = zeros(1,ntime+1);
ii = 0;
for ti = 0:deltt:ntime*deltt
	ii=ii+1;
	for i=1:2:100
		esol(ii) = esol(ii)+(1/i)*exp(-i*i*pi*pi*ti/25)*sin(i*pi/2);
	end 
end 
esol=100-(100*4/pi)*esol;

%-----------------------------
% plot the solution at nodes 8 
%-----------------------------------

time = 0:deltt:ntime*deltt;
plot(time,sol(1,:),'*',time,esol,'-');
xlabel('Time')
ylabel('Solution at nodes')
