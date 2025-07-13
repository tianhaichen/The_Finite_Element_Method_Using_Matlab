%-------------------------------------------------------
% Example 7.5.1
% to solve natural frequency of 1-d bar structure
%
% Problem description                                                        
%   Find the natural frequency of a bar structure           
%   as shown in Fig.  7.5.1          
%                                                                            
% Variable descriptions                                                      
%   k = element stiffness matrix                                             
%   kk = system stiffness matrix                                             
%   ff = system force vector                                                 
%   index = a vector containing system dofs associated with each element     
%   gcoord = global coordinate matrix
%   prop = element property matrix
%   nodes = nodal connectivity matrix for each element
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------%

%--------------------
%  control input data
%--------------------

clc
clear
nel = 4;            % number of elements
nnel = 2;           % number of nodes per element 
ndof = 1;           % number of dofs per node
nnode = 5;          % total number of nodes in system
sdof = nnode*ndof;  % total system dofs

%---------------------------------
%  nodal coordinates 
%---------------------------------

gcoord = [0,0;1,0;2,0;3,0;4,0];   % x coordinates of nodes

%----------------------------------
% material and geometric properties
%----------------------------------

prop= [200e9,0.001,7860];      % elastic modulus and cross-section density


%--------------------------------
% nodal connectivity
%-------------------------------

nodes(1,1)=1; nodes(1,2)=2;     % nodes associated with element 1
nodes(2,1)=2; nodes(2,2)=3;
nodes(3,1)=3; nodes(3,2)=4;
nodes(4,1)=4; nodes(4,2)=5;

%-------------------------
% applied constraints
%-------------------------

% bcdof=[1];    %constrainted
% bcval=[0];    % described value is 0

%---------------------------------
% initialization to zero
%---------------------------

kk = zeros(sdof,sdof);      % system stiffness matrix
mm = zeros(sdof,sdof);      % system mass matrix
index = zeros(nnel*ndof,1);     % index vector 

%-----------------
% loop for elements 
%------------------

for iel=1:nel   % loop for the total number of elements
    nd(1)=nodes(iel,1);     % 1st connected node for the (iel)-th element
    nd(2)=nodes(iel,2);

    x1=gcoord(nd(1),1); 
    x2=gcoord(nd(2),1); 

    leng = (x2-x1);  % element length

    el = prop(1);     % extract elastic modulus
    area = prop(2);   % extract cross-sectional area 
    rho = prop(3);      % extract mass dencity
    index = feeldof(nd,nnel,ndof); % extract system dofs for the element 
    ipt = 1;
    [k,m]=fetruss1(el,leng,area,rho,ipt);   
    kk=feasmbl1(kk,k,index);        % assemble into system matrix
    mm=feasmbl1(mm,m,index);
end

%------------------------------------
% solve for eigenvalues
%------------------------------------

fsol=eig(kk,mm);
fsol=sqrt(fsol);

%----------------------------------
% print fem solutions
%------------------------

num=1:1:sdof;
displ = [num' fsol]

