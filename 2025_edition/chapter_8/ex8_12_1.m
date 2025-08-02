%-----------------------------------------------------------------------
% Example 8.12.1
% to solve a one end fixed beam natural frequency and eigenvector using 
% Hermitian beam elements
%
% Problem description
%    Find the natural frequency and eigenvector of a one end fixed beam whose length is
%    1.27m. The beam has also elastic modulus of 10^7Pa and
%    moment of inertia of cross-section xi = 1/12.
%    Use 2 elements for one half of the beam due to symmetry.
%
% Variable descriptions
%   k = element stiffness matrix
%   m = element mass matrix
%   kk = system stiffness matrix
%   mm = system mass matrix
%   ff = system force vector
%   index = a vector containing system dofs associated with each element
%   bcdof = a vector containing dofs associted with boundary conditions
%   bcval = a vector containing boundary condition values associated with
%           the dofs in 'bcdof'
%----------------------------------------------------------------------------

clc
clear
nel = 2;                % number of elements
nnel = 2;               % number of nodes per element
ndof = 2;               % number of dofs per node
nnode = (nnel-1)*nel+1; % total number of nodes in system
sdof = nnode*ndof;      % total system dofs

el = 100;              % elastic modulus
xi = 0.01;              % moment of inertia of cross-section
tleng = 1.0;             % length of the beam
leng = tleng/nel;          % element length of equal size
area = 1;               % cross-sectional area of the beam
rho = 1;                % mass density 
ipt = 1;                % option for mass matrix (not used for static analysis)

bcdof = [1 2];             % first,12th dof is constrained
bcval = [0 0];              % value is 0

ff = zeros(sdof,1);     % initialization of system force vector
mm = zeros(sdof,sdof);     % initialization of system mass matrix
kk = zeros(sdof,sdof);  % initialization of system matrix 
index = zeros(nnel*ndof,1); %initialization of index vector

for iel = 1:nel
    index = feeldof1(iel,nnel,ndof);    % extract system dofs associated with element
    [k,m] = febeam1(el,xi,leng,area,rho,ipt);   % compute element stiffness matrix
    kk = feasmbl1(kk,k,index);
    mm = feasmbl1(mm,m,index);
end

% [kk,ff] = feaplyc2(kk,ff,bcdof,bcval); % apply the boundary conditions
% [mm,ff] = feaplyc2(kk,ff,bcdof,bcval); % apply the boundary conditions
[Omega,Phi,ModF] = femodal(mm,kk,ff)

% print both exact and fem solutions
