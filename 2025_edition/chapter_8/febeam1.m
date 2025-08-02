function [k,m] = febeam1(el,xi,leng,area,rho,ipt)
% -----------------------------------------------------------------------
% Purpose:
%   Stiffness and mass matrices for Hermitian beam element node dof
%       (v_1 theta_1 v_2 theta_2)
%
% Synopsis:
%   [k,m] = = febeam1(el,xi,leng,area,rho,itp)
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
c = el*xi/leng^3;
k11 =12;
k12 = 6*leng;
k22 = 4*leng^2;
k24 = 2*leng^2;
k=c*[k11 k12 -k11 k12;
    k12 k22 -k12 k24;
    -k11 -k12 k11 -k12;
    k12 k24 -k12 k22];

% consistent mass matrix 
if ipt ==1
mm=rho*area*leng/420;
m11 = 156;
m12 = 22*leng;
m13 = 54;
m14 = -13*leng;
m22 = 4*leng^2;
m24 = -3*leng^2;
m = mm*[m11 m12 m13 m14;
        m12 m22 -m14 m24;
        m13 -m14 m11 -m12;
        m14 m24 -m12 m22];
%
% lumper mass matrix
elseif ipt ==2
    m=zeros(4,4);
    mass=rho*area*leng;
    m=diag(mass/2*[1 0 1 0])
%
% diagonal mass matrix
else
    m=zeros(4,4);
    mass=rho*area*leng/78;
    m=diag(mass*[39 leng^2 39 leng^2]);
end
end