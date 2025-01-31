function J = PhiBarJacobian(zbar,TinOd,sigmaSquared)
%PhiBarJacobian Return the Jacobian J(z̄) of the balancing transformation x = ̅Φ(z̄) = Φ(𝝋(z̄)) (at one point).
%
%   Usage:  J = PhiBarJacobian(zbar,TinOd,sigmaSquared)
%
%   Inputs:     zbar  - point at which to evaluate the transformation
%               TinOd - cell array containing input-normal/output-diagonal
%                       transformation coefficients
%        sigmaSquared - the coefficients of the squared singular value functions
%
%   Outputs:       J  - the value of the Jacobian of the balancing transformation
%
%   Description: The nonlinear balancing transformation is x = ̅Φ(z̄) =
%   Φ(𝝋(z̄), which composition of the input-normal/output-diagonal
%   transformation
%              x = Φ(z)
%                = T₁z + T₂(z⊗z) + ... + Td(z...⊗z)
%   and the scaling transformation
%              z = 𝝋(z̄) , where zᵢ = 𝜑ᵢ(z̄) = z̄ᵢ / √ ̅σᵢ(z̄ᵢ)
%
%   The Jacobian J(z̄) is the matrix giving the partial derivatives ∂
%   ̅Φ(z̄)/∂z̄ of the transformation ̅Φ(z̄). Since this transformation is
%   the composition of two transformations, the Jacobian can be computed
%   using the chain rule as
%           J(z̄) = ∂ ̅Φ(z̄)/∂z̄
%                = [∂Φ(z)/∂z] [∂𝝋(z̄)/∂z̄]
%   where [∂Φ(z)/∂z] is the Jacobian of the input-normal/output-diagonal
%   transformation evaluated at the point z = 𝝋(z̄) and [∂𝝋(z̄)/∂z̄] is
%   the Jacobian of the scaling transformation evaluated at z̄. So in this
%   function, two main things need to be computed:
%       1) Given z̄, we need to compute z. We can do this via Newton
%          iteration, the same as in the function PhiBar().
%       2) We need the Jacobian [∂𝝋(z̄)/∂z̄] .
%   The function 𝝋(z̄) is given by
%           𝝋ᵢ(̄zᵢ) = z̄ᵢ / √ ̅σᵢ(z̄ᵢ)       ( = zᵢ )
%   and its Jacobian is the diagonal matrix
%           [∂𝝋(z)/∂z]ᵢⱼ = 
%            ( √ ̅σᵢ(z̄ᵢ) - z̄ᵢ ̅σᵢ'(z̄ᵢ)/(2√ ̅σᵢ(z̄ᵢ)) ) / ̅σᵢ(z̄ᵢ) for i=j, 0 else
%   However, we can't compute  ̅σᵢ(z̄ᵢ), only σᵢ(zᵢ). No problem; instead
%   of computing the diagonal matrix [∂𝝋(z)/∂z], we will compute its
%   inverse and then invert it, because the inverse of the Jacobian the
%   same as the Jacobian of the inverse! Fortunately, we already know how
%   to compute [∂𝝋⁻¹(z)/∂z], which was done in the function PhiBar(). Then
%   we simply use the fact that [∂𝝋(z)/∂z] = [∂𝝋⁻¹(z)/∂z]⁻¹. In fact,
%   [∂𝝋⁻¹(z)/∂z] will already be defined in the process of computing z
%   given z̄, and since it is a diagonal matrix, inversion is not expensive.
% 
%   References: [1]
%
%   Part of the NLbalancing repository.
%%

n = length(zbar);
dsigmaSquared = zeros(size(sigmaSquared) - [0 1]);
sigmaSquared = flip(sigmaSquared,2);

for i=1:n
    dsigmaSquared(i,:) = polyder(sigmaSquared(i,:));
end

%% Compute z given z̄ via Newton iteration
% Define functions for σ(z) and its derivative σ'(z); sigma can be defined
% as an anonymous function, but dsigma requires evaluating sigma(z) and
% indexing the components, which can't be done in one line as an anonymous
% function so I have to do it as a nested function
sigma = @(z) arrayfun(@(i) real(polyval(sigmaSquared(i,:), z(i))^(1/2)), 1:n).';
    function ds = dsigma(z)
        s = sigma(z);
        ds = arrayfun(@(i) polyval(sigmaSquared(i,:), z(i)) / (2*s(i)), 1:n).';
    end

% Define function and Jacobian for Newton iteration
varphiInv = @(z) z .* sqrt(sigma(z)); % 𝝋⁻¹(z)
dvarphiInv = @(z) diag(sqrt(sigma(z)) + z .* dsigma(z) ./ (2 * sqrt(sigma(z)))); % [∂𝝋⁻¹(z)/∂z]

% Solve for z using Newton iteration
z = newtonIteration(zbar, varphiInv, dvarphiInv);

%% Evaluate J(z̄) = ∂ ̅Φ(z̄)/∂̄z: composition of Jacobians [∂Φ(z)/∂z]*[∂𝝋(z̄)/∂z̄]
% The key is that the Jacobian of the inverse is the inverse of the Jacobian,
% so [∂𝝋(z̄)/∂z̄] = [∂𝝋⁻¹(z)/∂z]⁻¹. So instead of multiplying, we divide
% by the dvarphiInv function which we have already defined!
J = jcbn(TinOd, z)  / dvarphiInv(z);

end