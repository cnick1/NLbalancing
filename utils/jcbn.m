function J = jcbn(T, z)
%jcbn Return the Jacobian of the transformation given by the T coefficients at z.
%
%   Usage:  J = jcbn(T, z)
%
%   Inputs: T - cell array containing transformation coefficients
%           z - point at which to evaluate the Jacobian
%
%   Description: The Jacobian is the matrix giving the partial derivatives
%   ∂Φ(z)/∂z of the transformation Φ(z). We approximate the transformation
%   as a polynomial, so
% 
%              x = Φ(z)
%                = T{1}*z + T{2}*(z⊗z) + ... + T{d}*(z...⊗z)
% 
%   where the rows of the transformation coefficients are symmetrized. Then,
%   the Jacobian is given by
% 
%       J(z) = ∂Φ(z)/∂z = T{1} + 2*T{2}*(I⊗z) + ... + d*T{d}*(I...⊗z)
% 
%   Evaluating this explicitly is expensive, so this function uses the
%   Kronecker-vec identity to do this recursively and in a more efficient
%   manner.
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%

n = size(z, 1);
if isempty(T{1})
    J = 0;
else
    J = T{1};
end

zkm1 = 1;
for k = 2:length(T)
    % Using kron-vec identity (level-3 vs level-2 BLAS)
    zkm1 = kron(zkm1, z);
    % Need to iterate over n rows to apply kron-vec identity
    for j = 1:n
        J(j,:) = J(j,:) + k * zkm1.' * reshape(T{k}(j,:),n^(k-1),[]); %#ok<AGROW>
    end
end
end