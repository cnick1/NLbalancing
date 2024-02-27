function v = approxPastEnergy(f, g, h, eta, degree, verbose)
%approxPastEnergy  Compute the past energy function for a polynomial control-affine dynamical system.
%
%   Usage: v = approxPastEnergy(f,g,h,eta,d,verbose)
%
%   Inputs:
%       f,g,h   - cell arrays containing the polynomial coefficients
%                 for the drift, input, and output.
%                   • f must contain at least linear and quadratic coefficients
%                   • g must contain at least a linear input (B matrix)
%                   • h must contain at least a linear input (C matrix)
%       eta     - η=1-1/γ^2, where γ is the H∞ gain parameter. For open-loop
%                 balancing, use eta=0. For closed-loop (HJB) balancing, use
%                 eta=1. Any other value between -1 and ∞ corresponds to
%                 H∞ balancing.
%       degree  - desired degree of the computed energy function. A degree d
%                 energy function uses information from f,g,h up-to degree d-1.
%                 The default choice of d is lf+1, where lf is the degree of
%                 the drift.
%       verbose - optional argument to print runtime information
%
%   Output:
%       v       - cell array containing the polynomial energy function coefficients
%
%   Background: Computes a degree d polynomial approximation to the past energy function
%
%          E^-(x) = 1/2 ( v{2}'*(x⊗x) + ... + v{d}'*(...⊗x) )
%
%   for the polynomial control-affine system
%
%    \dot{x} = Ax + F2*(x⊗x) + F3*(x⊗x⊗x) + ...
%              + Bu + G1*(x⊗u) + G2*(x⊗x⊗u) + ...
%          y = Cx + H2*(x⊗x) + H3*(x⊗x⊗x) + ...
%
%   where eta = η=1-1/γ^2, where γ is the H∞ gain parameter. v{2} = vec(V2) = V2(:)
%   solves the Algebraic Riccati Equation
%
%    A'*V2 + V2*A + V2*B*B'*V2 - eta*C'*C = 0.
%
%   and the remaining v{i} solve linear systems arising from the Past H∞
%   Hamilton-Jacobi-Bellman Partial Differential Equation.
%
%   Details are in Section III.B of reference [1] or III.A of reference [2].
%
%   Requires the following functions from the KroneckerTools repository
%      KroneckerSumSolver
%      kronMonomialSymmetrize
%      LyapProduct
%
%   Authors: Jeff Borggaard, Virginia Tech
%            Nick Corbin, UCSD
%
%   License: MIT
%
%   Reference: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki, “Nonlinear
%               balanced truncation: Part 1—computing energy functions,” arXiv,
%               Dec. 2022. doi: 10.48550/ARXIV.2209.07645
%              [2] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” 2023.
%
%             See Algorithm 1 in [1].
%
%  Part of the NLbalancing repository.
%%

if (nargin < 6)
    verbose = false;
    if (nargin < 5)
        degree = length(f{1});
    end
end

% Print what type of energy function is being computed
if eta == 0
    message = sprintf('Computing open-loop balancing controllability energy function (η=%g ↔ γ=%g)', eta, 1 / sqrt(1 - eta));
    q = 0;
elseif eta == 1
    message = sprintf('Computing closed-loop balancing past energy function (η=%g ↔ γ=%g)', eta, 1 / sqrt(1 - eta));
    q = cellfun(@(x) x * (-1), h2q(h), 'un', 0);
else
    message = sprintf('Computing 𝓗∞ balancing past energy function (η=%g ↔ γ=%g)', eta, 1 / sqrt(1 - eta));
    q = cellfun(@(x) x * (-eta), h2q(h), 'un', 0);
end
if verbose
    disp(message)
end

% Rewritten by N Corbin to use ppr()
[v] = ppr(f, g, q, -1, degree, true, verbose);

end
