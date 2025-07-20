%-----------------------------------------------------------------------
% Example 8.10.2
% to find the natural frequencies and mode shapes of a using beam elements
% displacement degrees of freedom only
%
% Problem description
%   Find the natural frequencies of a clamped beam whose length is 1 m.      %
%   long. The beam has the cross-section of 0.02 m by 0.02 m and the mass density            %
%   is 1000 Kg/m^3. The elastic and shear modulus is 100 GPAa and         %
%   40 GPa, respectively. Use 4 elements. 
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
nel = 4;                % number of elements
nnel = 2;               % number of nodes per element
ndof = 3;               % number of dofs per node
nnode = (nnel-1)*nel+1; % total number of nodes in system
sdof = nnode*ndof;      % total system dofs

el = 100*10^9;              % elastic modulus
sh=40*10^9;      % shear modulus
tleng=1;         % total beam length
leng=tleng/nel;  % same size of beam elements
heig=0.02;       % height (or thickness) of the beam
width=0.02;      % width of the beam
rho=1000;        % mass density of the beam

bcdof=[1 2 3];

mm = zeros(sdof,sdof);     % initialization of system mass matrix
kk = zeros(sdof,sdof);  % initialization of system matrix 
index = zeros(nel*ndof,1); %initialization of index vector

for iel = 1:nel
    index = feeldof1(iel,nnel,ndof);    % extract system dofs associated with element
    [k,m] = febeam3(el,sh,leng,heig,width,rho);   % compute element stiffness matrix
    kk = feasmbl1(kk,k,index);
    mm = feasmbl1(mm,m,index);
end

[kn,mn] = feaplycs(kk,mm,bcdof);

fsol = eig(kn,mn);
fsol = sqrt(fsol)