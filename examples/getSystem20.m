function [f,g,h] = getSystem20(u,alpha,Nv,Nx,L,nu)
%getSystem20 Vlasov discretized with Hermite/Fourier
%
%   Usage:  [f,g,h] = getSystem20(u,alpha,Nv,Nx,L,nu)
%
%   Background: 
%
%   References: [1] 
%
%   Part of the PPR repository.
%%

% Define Parameters
NxTotal = 2*Nx+1;

%% Define A
% Define D
D = diag((2 * pi * 1i * (-Nx:Nx)) / L); D = sparse(D);

A_diag = kron(speye(Nv),D); 
A_offDiag = kron(diag(sqrt((1:Nv-1)/2),1),D) + kron(diag(sqrt((1:Nv-1)/2),-1),D);

n = 0:Nv-1; 
A_col = kron(diag((n .* (n - 1) .* (n - 2)) / ((Nv - 1) * (Nv - 2) * (Nv - 3))),speye(2*Nx+1));

A = -(u*A_diag + alpha*A_offDiag + nu * A_col);

%%% Define g(x)
%% Define B1 
Bstar = sqrt(2)/alpha * kron(diag(sqrt((1:Nv-1)),-1),speye(2*Nx+1));
W = dftmtx(2*NxTotal-1); 

phi1 = zeros(2*NxTotal-1,NxTotal); 
phi1(1:NxTotal,1:NxTotal) = eye(NxTotal); phi1 = sparse(phi1);

phi2 = zeros(NxTotal,2*NxTotal-1); 
phi2(1:NxTotal,(1:NxTotal)+Nx) = eye(NxTotal); phi2 = sparse(phi2);

B1 = Bstar *  kron(speye(Nv), 1/(2*NxTotal-1) * phi2 * W'*kr((W*phi1).',(W*phi1).').' * kron(speye(NxTotal),D));

%% Define B2
B2 = sparse(Nv*(2*Nx+1),2*Nx+1);
B2((2*Nx+2):2*(2*Nx+1),:) = 1/alpha^2 * D;

%% Define "K" -> like y=Cx output
DinvDiag = 1./((2 * pi * 1i * (1:Nx)) / L); 
DinvDiag = [-flip(DinvDiag),0,DinvDiag]; 
Dinv = diag(DinvDiag); Dinv = sparse(Dinv);

K = sparse(2*Nx+1,Nv*(2*Nx+1));
K(:,1:NxTotal) = alpha * Dinv.^2;

f = {A};
g = {B2,B1}; 
h = {K};


end

function X = kr(U,varargin)
%KR Khatri-Rao product.
%   kr(A,B) returns the Khatri-Rao product of two matrices A and B, of 
%   dimensions I-by-K and J-by-K respectively. The result is an I*J-by-K
%   matrix formed by the matching columnwise Kronecker products, i.e.
%   the k-th column of the Khatri-Rao product is defined as
%   kron(A(:,k),B(:,k)).
%
%   kr(A,B,C,...) and kr({A B C ...}) compute a string of Khatri-Rao 
%   products A o B o C o ..., where o denotes the Khatri-Rao product.
%
%   See also kron.
%   Version: 21/10/10
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
if ~iscell(U), U = [U varargin]; end
K = size(U{1},2);
if any(cellfun('size',U,2)-K)
    error('kr:ColumnMismatch', ...
          'Input matrices must have the same number of columns.');
end
J = size(U{end},1);
X = reshape(U{end},[J 1 K]);
for n = length(U)-1:-1:1
    I = size(U{n},1);
    A = reshape(U{n},[1 I K]);
    X = reshape(bsxfun(@times,A,X),[I*J 1 K]);
    J = I*J;
end
X = reshape(X,[size(X,1) K]);
end

