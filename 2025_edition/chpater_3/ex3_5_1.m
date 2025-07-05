%----------------------------------------------------------------------------
% EX3.5.1                                                              
% to solve the ordinary differential equation given as 
% au'' + bu' + cu = 1 ,0<x<1
% u(0)=0 and u(1)=1
% using 5 linear elements
%
% Variable descriptions
% k = element matrix
% f = element vector
% kk = system matrix
% index = a vector containing system dofs associated with each element
% bcdof = a vector containing dofs associated with boundary conditions
% bcval = a vector constaining boundary conditions values associated with 
%       the dofs in 'bcdof'
%----------------------------------------------------------------------------            

%-----------------------------------
% input data for control parameters
%-----------------------------------

clear
nnode = 50;          % total number of nodes in system
nel = nnode-1;      % number of elements
nnel = 2;           % number of nodes per element
ndof = 1;           % number of dofs per node

sdof = nnode*ndof;  % total system dofs

%-----------------------------------
% input data for nodal coordinate values
%-----------------------------------

% gcoord =[0.0 0.2 0.4 0.6 0.8 1.0];
gcoord=0:1/(nnode-1):1;
%---------------------------------------------------
% input data for nodal connectivity for each element
%---------------------------------------------------
 
% nodes=[1 2;        2 3;        3 4;        4 5;        5 6];
for i=1:nel
   for j=1:2
      nodes(i,j) = i+j-1;
   end
end

%---------------------------------------------
% input data for coefficients of the ODE
%---------------------------------------

acoef = 1;          % coefficient 'a' of the diff eqn
bcoef = -3;         % coefficient 'b' of the diff eqn 
ccoef = 2;          % coefficient 'c' of the diff eqn

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof(1) = 1;           % first node is constrained
bcval(1) = 0;           % whose described value is 0
bcdof(2) = nnode;           % 6th node is constrained
bcval(2) = 0;           % whose described value is 0

%----------------------------------------
% initialization of matrices and vectors
%----------------------------------------

ff = zeros(sdof,1);         % initialization of system force vector
kk = zeros(sdof,sdof);      % initialization of system matrix
index = zeros(nnel*ndof,1); % initialization of index vector

%----------------------------------------------------------------
% computation of element matrices and vectors and their assembly
%----------------------------------------------------------------

for iel = 1:nel         % loop for the total number of elements
    nl = nodes(iel,1); nr = nodes(iel,2); % extract nodes for iel-th element
    xl = gcoord(nl); xr = gcoord(nr);       % extract nodal coord values for the element
    eleng = xr - xl;                        % element length
    index = feeldof1(iel,nnel,ndof);        % extract system dofs associated with element
    k = feode2l(acoef,bcoef,ccoef,eleng);		% compute element matrix 
	f = fef1l(xl,xr);
    [kk,ff] = feasmbl2(kk,ff,k,f,index);    % assemble element matrices and vectors
end
%
%--------------------------
% apply boundary conditions
%--------------------------
[kk,ff] = feaplyc2(kk,ff,bcdof,bcval);
%
%-----------------------------
% solve the matrix equation
%-----------------------------
fsol=kk\ff;
%
%--------------------------
% analytical solution 
%--------------------------
c1=0.5/exp(1);
c2=-0.5*(1+1/exp(1));
for i=1:nnode
   x=gcoord(i);
   esol(i)=c1*exp(2*x)+c2*exp(x)+1/2;
end
%
%-------------------------------
% print both exact and fem solutions
%------------------------------
num=1:1:sdof;
retults=[num' fsol esol']
plot(num',fsol,'*',num',esol','-');
xlabel('num');
ylabel('solution');