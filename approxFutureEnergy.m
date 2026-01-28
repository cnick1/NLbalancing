function [w, K] = approxFutureEnergy(f, g, h, nvp)
%approxFutureEnergy  Compute the future energy function for a polynomial control-affine dynamical system.
%
%   Usage: w = approxFutureEnergy(f,g,h)
%
%   Inputs:
%       f,g,h   - cell arrays containing the polynomial coefficients
%                 for the drift, input, and output.
%                   • f must contain at least a linear drift  (A matrix)
%                   • g must contain at least a linear input  (B matrix)
%                   • h must contain at least a linear output (C matrix)
%
%     Optional name/value pair inputs:
%           eta - η=1-1/γ², where γ is the H∞ gain parameter.
%                  • For open-loop balancing, use eta=0.
%                  • For closed-loop (HJB) balancing, use eta=1.
%                 Any other value between -1 and ∞ corresponds to H∞
%                 balancing. The default is 0 for open-loop balancing.
%        degree - desired degree of the computed energy function. A degree d
%                 energy function uses information from f,g,h up-to degree d-1.
%                 The default choice of d is lf+1, where lf is the degree of
%                 the drift.
%       verbose - optional argument to print runtime information
%
%   Output:
%       w       - cell array containing the polynomial energy function coefficients
%       K       - the gain coefficients corresponding to the optimal controller
%                 given by the future energy function (optional)
%
%   Description: For control-affine dynamics ẋ = f(x) + g(x) u, y = h(x) we
%   seek an approximation of the H∞ future energy function
%
%           E⁺(x) = minᵤ J(x,u) = ½∫ ||y||² + ||u||²/η dt
%
%   where η=1-1/γ² and γ is the H∞ gain parameter. The solution is given by the
%   solution to the HJB PDE
%
%           0 = 𝜕ᵀE⁺(x)/𝜕x f(x) - η/2 𝜕ᵀE⁺(x)/𝜕x g(x) gᵀ(x) 𝜕E⁺(x)/𝜕x + ½ h(x)ᵀ h(x)
%
%   A local approximation can be computed using the method of Al'brekht [2],
%   i.e. we compute the Taylor expansions:
%
%           E⁺(x) = 1/2 ( w₂ᵀ(x ⊗ x) + w₃ᵀ(x ⊗ x ⊗ x) + ... +   wᵈᵀ(... ⊗ x) )
%
%   based on the Taylor expansions for the dynamics, written as
%
%           ẋ = A x + F₂ (x ⊗ x) + F₃ (x ⊗ x ⊗ x) + ...
%               + B u + G₁ (x ⊗ u) + G₂ (x ⊗ x ⊗ u) + ...
%           y = C x + H₂ (x ⊗ x) + H₃ (x ⊗ x ⊗ x)
%
%   Inserting all these polynomial expressions into the HJB PDEs (1) and (2)
%   leads to equations for the energy function coefficients w₂, w₃,..., wᵈ. The
%   first coefficient w₂ = vec(W₂) = W₂(:) solves the H∞ Algebraic Riccati
%   Equation
%
%           Aᵀ W₂ + W₂ A + Cᵀ C - η W₂ B Bᵀ W₂ = 0,
%
%   The remaining wᵢ solve linear systems arising from (1).
%
%   Details are in Section III.B of reference [1] or III.A of reference [2], and
%   in Algorithm 1 in both references.
%
%   The solution is computed using the ppr() function in the PPR repository and
%   requires the following functions from the KroneckerTools repository:
%      KroneckerSumSolver
%      kronMonomialSymmetrize
%      LyapProduct
%      h2q
%
%   Authors: Rewritten by Nick Corbin, UCSD to use the PPR package [3]
%            Original version by Jeff Borggaard, Virginia Tech
%
%   License: MIT
%
%   References: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki,
%               “Scalable computation of energy functions for nonlinear balanced
%               truncation,” Computer Methods in Applied Mechanics and
%               Engineering, vol. 427, p. 117011, Jul. 2024, doi:
%               10.1016/j.cma.2024.117011
%               [2] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗∞
%               energy functions for polynomial control-affine systems," IEEE
%               Transactions on Automatic Control, pp. 1–13, 2024, doi:
%               10.1109/tac.2024.3494472
%               [3] N. A. Corbin and B. Kramer, “Computing solutions to the
%               polynomial-polynomial regulator problem,” in 2024 IEEE 63rd
%               Conference on Decision and Control (CDC), IEEE, Dec. 2024, pp.
%               2689–2696. doi: 10.1109/cdc56724.2024.10885897.
%
%   Part of the NLbalancing repository.
%
%   See also: ppr
%%
arguments
    f
    g
    h
    nvp.eta = 0
    nvp.degree = length(f)
    nvp.verbose = false
    nvp.r = size(f{1},1)
end

% Print what type of energy function is being computed
if nvp.eta == 0
    message = sprintf('Computing open-loop balancing observability energy function (η=%g ↔ γ=%g)', nvp.eta, 1 / sqrt(1 - nvp.eta));
    nvp.eta = Inf; % Need 1/eta to be zero, and if eta = 0 this doesn't work. Basically R = R^-1 = 0 is what we need in ppr()
elseif nvp.eta == 1
    message = sprintf('Computing closed-loop balancing future energy function (η=%g ↔ γ=%g)', nvp.eta, 1 / sqrt(1 - nvp.eta));
else
    message = sprintf('Computing 𝓗∞ balancing future energy function (η=%g ↔ γ=%g)', nvp.eta, 1 / sqrt(1 - nvp.eta));
end
if nvp.verbose
    disp(message)
end

options.verbose = nvp.verbose; options.reducedDimension = nvp.r;
[w, K] = ppr(f, g, h2q(h), 1/nvp.eta, nvp.degree, options);

end
