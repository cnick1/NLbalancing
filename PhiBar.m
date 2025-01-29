function x = PhiBar(zbar,TinOd,sigmaSquared)
%PhiBar Return the balancing transformation x = ̅Φ(z̄) = Φ(𝝋(z̄)) (at one point).
%
%   Usage:  x = PhiBar(zbar,TinOd,sigmaSquared)
%
%   Inputs:     zbar  - point at which to evaluate the transformation
%               TinOd - cell array containing input-normal/output-diagonal
%                       transformation coefficients
%        sigmaSquared - the coefficients of the squared singular value functions
%
%   Outputs:       x  - the value of x = ̅Φ(z̄) = Φ(𝝋(z̄)
%
%   Description: The nonlinear balancing transformation is x = ̅Φ(z̄) =
%   Φ(𝝋(z̄), which composition of the input-normal/output-diagonal
%   transformation x = Φ(z) and the scaling transformation z = 𝝋(z̄).
%   The balancing transformation puts the energy functions in the form
%           𝓔⁻( ̅Φ(z̄)) = 1/2 z̄ᵀ Σ⁻¹(z̄) z̄
%           𝓔⁺( ̅Φ(z̄)) = 1/2 z̄ᵀ  Σ(z̄)  z̄
%   where Σ(z̄) is the diagonal matrix of singular value functions ̅σᵢ(z̄ᵢ).
%   There is basically one major thing to compute:
%       1) Given z̄, we need to compute z.
%   This requires solving the equations zᵢ √σᵢ(zᵢ) = z̄ᵢ for zᵢ. In other
%   words, while the transformation is z = 𝝋(z̄), what we really have
%   access to is the inverse transformation z̄ = 𝝋⁻¹(z), so we need to
%   solve z̄ = 𝝋⁻¹(z) for z given z̄. We can do this via Newton iteration.
%   The function 𝝋⁻¹(z) is given by
%           𝝋ᵢ⁻¹(zᵢ) = zᵢ √σᵢ(zᵢ)
%   and its Jacobian is the diagonal matrix
%           [∂𝝋⁻¹(z)/∂z]ᵢⱼ = √σᵢ(zᵢ) + zᵢ σᵢ'(zᵢ)/(2√σᵢ(zᵢ)) for i=j, 0 else
%   Once we get z, we can get x by evaluating the input-normal transformation
%   for the point z as
%           x = Φ(z) = T₁z + T₂(z⊗z) + ... + Td(z...⊗z)
%   using the kronPolyEval(TinOd,z) function.
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

%% Evaluate x = ̅Φ(z̄): compute x given z
x = kronPolyEval(TinOd, z);

end
