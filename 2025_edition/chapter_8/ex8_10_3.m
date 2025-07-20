%-----------------------------------------------------------------------
% Example 8.10.2
% to find the natural frequencies for a 2-d frame using frame elements       %
%                                                                            %
% Problem description                                                        %
%   Find the natural frequencies of a frame of L-shape which is made of      %
%   two beams of length of 1 m each. Both beams have                         %
%   cross-sections of 0.01 m by 0.01 m. The elastic modulus is 100 GPa.      %
%   The beam has mass density of 1000 Kg/m^3.  Use 10 elements. 
%
% Variable descriptions
%   k = element stiffness matrix
%   m = element mass matrix
%   kk = system stiffness matrix
%   mm = system mass matrix
%   index = a vector containing system dofs associated with each element
%   bcdof = a vector containing dofs associted with boundary conditions
%   bcval = a vector containing boundary condition values associated with
%           the dofs in 'bcdof'
%----------------------------------------------------------------------------
addpath(genpath('D:\The_Finite_Element_Method_Using_Matlab\2025_edition'));
clc
clear
nel = 10;                % number of elements
nnel = 2;               % number of nodes per element
ndof = 3;               % number of dofs per node
nnode = (nnel-1)*nel+1; % total number of nodes in system
sdof = nnode*ndof;      % total system dofs
coord=[0,0;
    0,0.2;
    0,0.4;
    0,0.6;
    0,0.8;
    0,1.0;
    0.2,1;
    0.4,1;
    0.6,1;
    0.8,1;
    1,1];
x=coord(:,1);
y=coord(:,2);

el = 100*10^9;              % elastic modulus
area=0.0001;        % cross-sectional area
xi=8.3333*10^(-10); % moment of inertia of cross-section
rho=1000;           % mass density per volume (dummy value for static analysis)
bcdof=[1 2 3];

mm = zeros(sdof,sdof);     % initialization of system mass matrix
kk = zeros(sdof,sdof);  % initialization of system matrix 
index = zeros(nel*ndof,1); %initialization of index vector

for iel = 1:nel
    index = feeldof1(iel,nnel,ndof);    % extract system dofs associated with element
    node1=iel;      % starting node number for element 'iel'
    node2=iel+1;    % ending node number for element 'iel'

    x1=x(node1); y1=y(node1); % x and y coordinate values of 'node1'
    x2=x(node2); y2=y(node2); % x and y coordinate values of 'node2'

    leng=sqrt((x2-x1)^2+(y2-y1)^2); % length of element 'iel'

    if (x2-x1)==0  % compute the angle between the local and global axes
        beta=pi/2;
    else
        beta=atan((y2-y1)/(x2-x1));
    end
    [k,m]=feframe2(el,xi,leng,area,rho,beta,1);   % compute element stiffness matrix
    kk = feasmbl1(kk,k,index);
    mm = feasmbl1(mm,m,index);
end

[kn,mn] = feaplycs(kk,mm,bcdof);

fsol = eig(kn,mn);
fsol = sqrt(fsol)