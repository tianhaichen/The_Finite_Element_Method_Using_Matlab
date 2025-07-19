function [k,m] = febeam4(el,xi,leng,sh,heig,rho,ipt)
% -----------------------------------------------------------------------
% Purpose:
%   Stiffness and mass matrices for mixed beam element bending moment
% and deflection as nodal degrees of freedom
% nodal dof   (M_1 v_1 M_2 v_2)
%
% Synopsis:
%   [k,m] = febeam4(el,xi,leng,sh,area,rho,ipt)
%
% Variable Description:
%   k = element stiffness matrix (size of 4x4)
%   m = element mass matrix (size of 4x4)
%   el = elastic modulus
%   xi = second moment of inertia of cross-section
%   leng = element length
%   heig - beam thickness
%   rho = mass density (mass per unit volume)
%   sh = shear modulus
%   ipt = 1: consistent mass matrix
%         2: lumped mass matrix
%         ohterwise: diagonal mass matrix
%________________________________________________________________________

% stiffness matrix
if sh==0
    % thin beam(no shear deformation)
    k11 = leng/(3*el*xi);
    k12 = 1/leng;
    k13 = leng/(6*el*xi);
    k24 = 0;
    k=[k11 k12 k13 -k12;
        k12 0 -k12 k24;
        k13 -k12 k11 k12;
        -k12 k24 k12 0];
else
    % thick beam (include shear deformation)
    a = 6/(5*sh*leng*heig);
    k11 = 1/(3*el*xi);
    k12 = 1/leng;
    k13 = 1/(6*el*xi);
    k=[k11+a k12 k13-a -k12;
        k12 0 -k12 0;
        k13-a -k12+a k11 k12;
        -k12 0 k12 0];
end



% consistent mass matrix
if ipt ==1
    % lumper mass matrix
    m=zeros(4,4);
    mass=rho*heig*leng/2;
    m=diag(mass*[0 1 0 1]);
    % diagonal mass matrix
else
    m=zeros(4,4);
    mass=rho*heig*leng/2;
    m=diag(mass*[1 1 1 1]);
end
end