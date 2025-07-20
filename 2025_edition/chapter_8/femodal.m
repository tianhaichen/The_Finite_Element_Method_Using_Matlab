function [Omega,Phi,ModF] = femodal(M,K,F)
%-------------------------------------------------------------------
%
%  Purpose:
%     The function subroutine femodal.m calculates modal parameters
%     for a given structural system. It calculates natural frequency
%       and eigenvector.The eigenvectors are normalized so that the 
%       modal mass matrix becomes an identity matrix.
%
%  Synopsis:
%   [Omega,Phi,ModF] = femodal(M,K,F)
%
%  Variables Description:
%    Input parameters - M,K -Mass and stiffness matrices
%                       F - Input or forcing influence matrix
%    Output parameters - Omega - Natural frequency(rad/sec) in
%                       ascending order
%                       -Phi - Modal matrix with each column corresponding
%                       to the eigenvector.
%                       ModF - Modal input patrices.
%----------------------------------------------------------------------
disp('  ')
disp('Please wait!! - The job is being performed.')

%---------------------------------------------------------------------
%  Solve the eigenvalue problem and normlized the eigenvectors
%---------------------------------------------------------------------
[n,n] = size(M);
[n,m] = size(F);
[V,D] =eig(K,M);
[lambda,k] = sort(diag(D));     % sort the eigenvalues and eigenvecotors
% in ascending order
V = V(:,k);
Factor = diag(V'*M*V);
Phi = V*inv(sqrt(diag(Factor)));    % Eigenvectors are normalized

Omega = diag(sqrt(Phi'*K*Phi));   % natural frequency in ascending order
ModF = Phi'*F;
end