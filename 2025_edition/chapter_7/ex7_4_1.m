%-------------------------------------------------------
% Example 7.4.1
% to solve static 2-d truss structure
%
% Problem description                                                        
%   Find the deflection and stress of the truss made of two members           
%   as shown in Fig. 7.4.1.           
%                                                                            
% Variable descriptions                                                      
%   k = element stiffness matrix                                             
%   kk = system stiffness matrix                                             
%   ff = system force vector                                                 
%   index = a vector containing system dofs associated with each element     
%   gcoord = global coordinate matrix
%   disp = nodal displacement vector
%   elforce = element force vector
%   eldisp = element nodal displacement
%   stress = stress vector for every element
%   elprop = element property matrix
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
nel = 2;            % number of elements
nnel = 2;           % number of nodes per element 
ndof = 2;           % number of dofs per node
nnode = 3;          % total number of nodes in system
sdof = nnode*ndof;  % total system dofs

%---------------------------------
%  nodal coordinates 
%---------------------------------

gcoord = [0,0;10.0,0;0.0,10];   % x,y coordinates of nodes

%----------------------------------
% material and geometric properties
%----------------------------------

elprop(1,1) = 300e5;    
elprop(1,2) = 0.4;      % elastic modulus and cross-section of element
elprop(2,1) = 30e5;
elprop(2,2) = 0.5;

%--------------------------------
% nodal connectivity
%-------------------------------

nodes(1,1)=1; nodes(1,2)=2;     % nodes associated with element 1
nodes(2,1)=2; nodes(2,2)=3;

%-------------------------
% applied constraints
%-------------------------

bcdof=[1,2,5,6];    %constrainted
bcval=[0,0,0,0];    % described value is 0

%---------------------------------
% initialization to zero
%---------------------------

ff = zeros(sdof,1);     % system force vector
kk = zeros(sdof,sdof);
index = zeros(nnel*ndof,1);     % index vector 
elforce = zeros(nnel*ndof,1);   % element force vector 
eldisp = zeros(nnel*ndof,1);    % element nodal displacement vector 
k = zeros(nnel*ndof,nnel*ndof); % element stiffness matrix 
stress = zeros(nel,1);          % stress vector for every element 

%--------------------
% applied nodal force 
%--------------------

ff(4) = -1000;      % 2nd node has 1000lb in downward direction 

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
    el = elprop(iel,1);     % extract elastic modulus
    area = elprop(iel,2);   % extract cross-sectional area     
    index = feeldof(nd,nnel,ndof); % extract system dofs for the element 
    k=fetruss2(el,leng,area,0,beta,1);
    kk=feasmbl1(kk,k,index);        % assemble into system matrix
end

%-----------------------------------------------------
% apply constraints and solve the matrix
%---------------------------------------

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval); % apply the boundary conditions

disp = kk\ff;    % solve the matrix equation to find nodal displacements

%--------------------------------------
% post computation for stress calculation
%----------------------------------------

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
    el = elprop(iel,1);     % extract elastic modulus
    area = elprop(iel,2);   % extract cross-sectional area     
    index = feeldof(nd,nnel,ndof);  % extract system dofs for the element 
    k=fetruss2(el,leng,area,0,beta,1);
    
    for i=1:(nnel*ndof)
        eldisp(i) = disp(index(i));
    end
    elforce = k*eldisp;     % element force vector
    stress(iel) = sqrt(elforce(1)^2+elforce(2)^2)/area;  % stress calculation
    
    if((x2-x1)*elforce(3))<0
        stress(iel)=-stress(iel);
    end
end

%----------------------------------
% print fem solutions
%------------------------

num=1:1:sdof;
displ = [num' disp]

numm=1:1:nel;
stresses = [numm' stress]