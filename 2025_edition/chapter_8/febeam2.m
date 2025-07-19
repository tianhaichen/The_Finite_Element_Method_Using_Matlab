function [k,m] = febeam2(el,xi,leng,sh,area,rho,ipt)
% -----------------------------------------------------------------------
% Purpose:
%   Stiffness and mass matrices for C^0 beam element 
% nodal dof   (v_1 theta_1 v_2 theta_2)
%
% Synopsis:
%   [k,m] = febeam2(el,xi,leng,sh,area,rho,ipt)
%
% Variable Description:
%   k = element stiffness matrix (size of 4x4)
%   m = element mass matrix (size of 4x4)
%   el = elastic modulus
%   xi = second moment of inertia of cross-section
%   leng = element length
%   area = area of beam cross-section
%   rho = mass density (mass per unit volume)
%   sh = shear modulus
%   ipt = 1: consistent mass matrix
%         2: lumped mass matrix
%         ohterwise: diagonal mass matrix
%________________________________________________________________________

% stiffness matrix
c = el*xi/leng;
d = (5/6)*sh*area/(4*leng);
k11 =4*d;
k12 = 2*d*leng;
k22 = c+d*leng^2;
k24 = -c+d*leng^2;
k=[k11 k12 -k11 k12;
    k12 k22 -k12 k24;
    -k11 -k12 k11 -k12;
    k12 k24 -k12 k22];

% consistent mass matrix 
if ipt ==1
mm=rho*area*leng/420;
m11 = 2;
m12 = 0;
m13 = 1;
m14 = 0;
m22 = 0;
m24 = 0;
m = mm*[m11 m12 m13 m14;
        m12 m22 -m14 m24;
        m13 -m14 m11 -m12;
        m14 m22 -m12 m22];
%
% lumper mass matrix
%
% diagonal mass matrix
else
    m=zeros(4,4);
    mass=rho*area*leng/2;
    m=diag(mass*[1 0 1 0]);
end
end