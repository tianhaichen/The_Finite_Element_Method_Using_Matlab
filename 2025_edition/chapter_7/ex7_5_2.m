%-------------------------------------------------------
% Example 7.5.2
% to solve natural frequency of 2-d truss structure
%
% Problem description                                                        
%   Find the natural frequency of a truss structure           
%   as shown in Fig. 7.4.2.           
%                                                                            
% Variable descriptions                                                      
%   k = element stiffness matrix                                             
%   m = element mass matrix
%   kk = system stiffness matrix                                             
%   mm = system mass vector                                                 
%   index = a vector containing system dofs associated with each element     
%   gcoord = global coordinate matrix
%   prop = element property matrix
%   nodes = nodal connectivity matrix for each element
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%---------------------------------------------------------------------%---%

%--------------------
%  control input data
%--------------------

clc
clear
nel = 9;            % number of elements
nnel = 2;           % number of nodes per element 
ndof = 2;           % number of dofs per node
nnode = 6;          % total number of nodes in system
sdof = nnode*ndof;  % total system dofs

%---------------------------------
%  nodal coordinates 
%---------------------------------

gcoord = [0,0;4.0,0;
    4.0,3;8.0,0;
    8.0,3.0;12.0,0];   % x,y coordinates of nodes

%----------------------------------
% material and geometric properties
%----------------------------------

prop(1,1) = 200e9;    
prop(1,2) = 0.0025;      % elastic modulus and cross-section of element
prop(1,3) = 7860;       %density

%--------------------------------
% nodal connectivity
%-------------------------------

nodes(1,1)=1; nodes(1,2)=2;     % nodes associated with element 1
nodes(2,1)=1; nodes(2,2)=3;
nodes(3,1)=2; nodes(3,2)=3;
nodes(4,1)=2; nodes(4,2)=4;
nodes(5,1)=3; nodes(5,2)=4;
nodes(6,1)=3; nodes(6,2)=5;
nodes(7,1)=4; nodes(7,2)=5;
nodes(8,1)=4; nodes(8,2)=6;
nodes(9,1)=5; nodes(9,2)=6;

%-------------------------
% applied constraints
%-------------------------

bcdof=[1,2,12];    %constrainted
bcval=[0,0,0];    % described value is 0

%---------------------------------
% initialization to zero
%---------------------------

kk = zeros(sdof,sdof);
mm = zeros(sdof,sdof);
index = zeros(nnel*ndof,1);     % index vector 

%-----------------
% loop for elements 
%------------------

for iel=1:nel   % loop for the total number of elements
    nd(1)=nodes(iel,1);     % 1st connected node for the (iel)-th element
    nd(2)=nodes(iel,2);

    x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);
    x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);

    leng = sqrt((x2-x1)^2+(y2-y1)^2);  % element length

    % angle between local and global axes
    if (x2-x1)==0
        if y2>y1
            beta=2*atan(1);
        else
            beta=-2*atan(1);
        end
    else
        beta = atan((y2-y1)/(x2-x1));
    end
    el = prop(1);     % extract elastic modulus
    area = prop(2);   % extract cross-sectional area   
    rho = prop(3);   
    index = feeldof(nd,nnel,ndof); % extract system dofs for the element 
    ipt=1;
    [k,m] = fetruss2(el,leng,area,rho,beta,ipt);

    kk=feasmbl1(kk,k,index);        % assemble into system matrix
    mm=feasmbl1(mm,m,index);
end

%-----------------------------------------------------
% apply constraints and solve the matrix
%---------------------------------------

[kk,mm]=feaplycs(kk,mm,bcdof); % apply the boundary conditions

fsol=eig(kk,mm);
fsol=sqrt(fsol);

%----------------------------
% print fem solutions
%----------------------------

num=1:1:sdof;
freqcy=[num' fsol]          % print natural frequency

%--------------------------------------------------------------------