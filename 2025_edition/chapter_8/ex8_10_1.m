%-----------------------------------------------------------------------
% Example 8.10.1
% to find the natural frequencies and mode shapes of a using Hermitian
%  beam elements
%
% Problem description
%    Find the natural frequencies and mode shapes of a free beam of length
%   1. It has a cross-section 1 by 1 and it has also mass dencity of 1.
%   The elastic modulus of the beam is 12.
%   Use 4 elements to model the whole beam such that nonsymmetric mode
%   shapes can be included. use also consistent mass matrices.
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
ndof = 2;               % number of dofs per node
nnode = (nnel-1)*nel+1; % total number of nodes in system
sdof = nnode*ndof;      % total system dofs

el = 12;              % elastic modulus
xi = 1/12;              % moment of inertia of cross-section
tleng = 1;             % length of a half of the beam
leng = tleng/nel;          % element length of equal size
area = 1;               % cross-sectional area of the beam
rho = 1;                % mass density 
ipt = 1;                % option for mass matrix (not used for static analysis)

mm = zeros(sdof,sdof);     % initialization of system mass matrix
kk = zeros(sdof,sdof);  % initialization of system matrix 
index = zeros(nnel*ndof,1); %initialization of index vector

for iel = 1:nel
    index = feeldof1(iel,nnel,ndof);    % extract system dofs associated with element
    [k,m] = febeam1(el,xi,leng,area,rho,ipt)   % compute element stiffness matrix
    kk = feasmbl1(kk,k,index);
    mm = feasmbl1(mm,m,index);
end

fsol = eig(kk,mm)
fsol = sqrt(fsol)