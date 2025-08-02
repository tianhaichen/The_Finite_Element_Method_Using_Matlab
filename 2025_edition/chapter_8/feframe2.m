function [k,m] = feframe2(el,xi,leng,area,rho,beta,ipt)
% -----------------------------------------------------------------------
% Purpose:
%   Stiffness and mass matrices for the 2-d frame element
%     nodal dof {u_1 v_1 theta_1 u_2 v_2 theta_2}
%
% Synopsis:
%   [k,m] = = feframe2(el,xi,leng,area,rho,beta,ipt)
%
% Variable Description:
%   k = element stiffness matrix (size of 4x4)
%   m = element mass matrix (size of 4x4)
%   el = elastic modulus
%   xi = second moment of inertia of cross-section
%   leng = element length
%   area = area of beam cross-section
%   rho = mass density (mass per unit volume)
%   ipt = 1: consistent mass matrix
%         2: lumped mass matrix
%         ohterwise: diagonal mass matrix
%________________________________________________________________________

% stiffness matrix
a=el*area/leng;
 c=el*xi/(leng^3);
 kl=[a   0           0            -a     0          0;...
     0   12*c        6*leng*c      0    -12* c      6*leng*c;...
     0   6*leng*c    4*leng^2*c    0    -6*leng*c   2*leng^2*c;...
    -a   0           0             a     0          0;...
     0  -12*c       -6*leng*c      0     12*c      -6*leng*c;...
     0   6*leng*c    2*leng^2*c    0    -6*leng*c   4*leng^2*c];

% rotation matrix

r=[ cos(beta)  sin(beta)  0   0          0         0;...
    -sin(beta)  cos(beta)  0   0          0         0;...
    0          0          1   0          0         0;...
    0          0          0   cos(beta)  sin(beta) 0;...
    0          0          0  -sin(beta)  cos(beta) 0;...
    0          0          0   0          0         1];
% stiffness matrix at global axis
k=r'*kl*r;

% consistent mass matrix

if ipt==1

    mm=rho*area*leng/420;
    ma=rho*area*leng/6;
    ml=[2*ma  0           0            ma    0            0;...
        0     156*mm      22*leng*mm   0     54*mm       -13*leng*mm;...
        0     22*leng*mm  4*leng^2*mm  0     13*leng*mm  -3*leng^2*mm;...
        ma    0           0            2*ma  0            0;...
        0     54*mm       13*leng*mm   0     156*mm      -22*leng*mm;...
        0    -13*leng*mm -3*leng^2*mm  0    -22*leng*mm   4*leng^2*mm];

    % lumped mass matrix

elseif ipt==2

    ml=zeros(6,6);
    mass=rho*area*leng;
    ml=mass*diag([0.5  0.5  0  0.5  0.5  0]);

    % diagonal mass matrix

else

    ml=zeros(6,6);
    mass=rho*area*leng;
    ml=mass*diag([0.5  0.5  leng^2/78  0.5  0.5  leng^2/78]);

end

% mass in the global system

m=r'*ml*r;
end