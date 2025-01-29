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
%       2) We need the Jacobian [∂𝝋(z̄)/∂z̄] *** Is the inverse of the Jacobian the same as the Jacobian of the inverse?
%   The function 𝝋(z̄) is given by
%           𝝋ᵢ(̄zᵢ) = z̄ᵢ / √ ̅σᵢ(z̄ᵢ)       ( = zᵢ )
%   and its Jacobian is the diagonal matrix
%           [∂𝝋(z)/∂z]ᵢⱼ = ( √ ̅σᵢ(z̄ᵢ) - z̄ᵢσᵢ'(z̄ᵢ)/(2√ ̅σᵢ(z̄ᵢ)) ) / ̅σᵢ(z̄ᵢ) for i=j, 0 else
%   We can leverage the relations ̅σᵢ(z̄ᵢ) = σᵢ(zᵢ) and zᵢ = z̄ᵢ / √ ̅σᵢ(z̄ᵢ)
%   to compute this as
%           [∂𝝋(z)/∂z]ᵢⱼ = ( √ σᵢ(zᵢ) - zᵢσᵢ'(zᵢ)/2 ) / σᵢ(zᵢ) for i=j, 0 else
%
%
%       zᵢ = z̄ᵢ / √ ̅σᵢ(z̄ᵢ)
%               or
%       z̄ᵢ = zᵢ √σᵢ(zᵢ)
%   So I need to solve z̄ᵢ = zᵢ √σᵢ(zᵢ) for zᵢ. In other words, given z̄ᵢ,
%   I need to solve the roots of the scalar equation
%       gᵢ(zᵢ) = zᵢ √σᵢ(zᵢ) - z̄ᵢ
%   or the vector equation
%       g(z) = z ⊙ √σ(z) - z̄


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

%% Compute z given z̄
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
f = @(z) z .* sqrt(sigma(z));
J = @(z) diag(sqrt(sigma(z)) + z .* dsigma(z) ./ (2 * sqrt(sigma(z))));

% Solve for z using Newton iteration
z = newtonIteration(zbar, f, J);

%% Compute Jacobian of scaling 𝝋(z̄), Jscal = ∂𝝋(z̄)/∂z̄
Jscal = diag(   ( sqrt(sigma(z)) - z .* dsigma(z)/2 )./sigma(z)   );

%% Evaluate J(z̄) = ∂ ̅Φ(z̄)/∂̄z: composition of Jacobians [∂Φ(z)/∂z]*[∂𝝋(z̄)/∂z̄]
J = jcbn(TinOd, z)  * Jscal;

end