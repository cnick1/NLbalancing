function w = approxFutureEnergy(f, g, h, eta, degree, verbose)
%approxFutureEnergy  Compute the future energy function for a polynomial control-affine dynamical system.
%
%   Usage: w = approxFutureEnergy(f,g,h,eta,d,verbose)
%
%   Inputs:
%       f,g,h   - cell arrays containing the polynomial coefficients
%                 for the drift, input, and output.
%                   • f must contain at least linear and quadratic coefficients
%                   • g must contain at least a linear input (B matrix)
%                   • h must contain at least a linear input (C matrix)
%       eta     - η=1-1/γ², where γ is the H∞ gain parameter. For open-loop
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
%       w       - cell array containing the polynomial energy function coefficients
%
%   Description: Computes a degree d polynomial approximation to the energy function
%
%          E^+(x) = 1/2 ( w{2}'*(x⊗x) + ... + w{d}'*(...⊗x) )
%
%   for the polynomial control-affine system
%
%    ẋ = Ax + F2*(x⊗x) + F3*(x⊗x⊗x) + ...
%              + Bu + G1*(x⊗u) + G2*(x⊗x⊗u) + ...
%          y = Cx + H2*(x⊗x) + H3*(x⊗x⊗x) + ...
%
%   where eta = η=1-1/γ², where γ is the H∞ gain parameter. w{2} = vec(W2) = W2(:)
%   solves the H∞ Algebraic Riccati Equation
%
%    A'*W2 + W2*A - eta*W2*B*B'*W2 + C'*C = 0,
%
%   and the remaining w{i} solve linear systems arising from the Future H∞
%   Hamilton-Jacobi-Bellman Partial Differential Equation.
%
%   Details are in Section III.B of reference [1] or III.A of reference [2].
%
%   Requires the following functions from the KroneckerTools repository:
%      KroneckerSumSolver
%      kronMonomialSymmetrize
%      LyapProduct
%      h2q
%
%   Authors: Jeff Borggaard, Virginia Tech
%            Nick Corbin, UCSD
%
%   License: MIT
%
%   Reference: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki,
%               “Scalable computation of energy functions for nonlinear
%               balanced truncation,” Computer Methods in Applied Mechanics
%               and Engineering, vol. 427, p. 117011, Jul. 2024, doi:
%               10.1016/j.cma.2024.117011
%              [2] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗∞
%               energy functions for polynomial control-affine systems,"
%               IEEE Transactions on Automatic Control, pp. 1–13, 2024,
%               doi: 10.1109/tac.2024.3494472
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
    message = sprintf('Computing open-loop balancing observability energy function (η=%g ↔ γ=%g)', eta, 1 / sqrt(1 - eta));
    eta = Inf; % Need 1/eta to be zero, and if eta = 0 this doesn't work. Basically R = R^-1 = 0 is what we need in ppr()
elseif eta == 1
    message = sprintf('Computing closed-loop balancing future energy function (η=%g ↔ γ=%g)', eta, 1 / sqrt(1 - eta));
else
    message = sprintf('Computing 𝓗∞ balancing future energy function (η=%g ↔ γ=%g)', eta, 1 / sqrt(1 - eta));
end
if verbose
    disp(message)
end

% Rewritten by N Corbin to use ppr()
options.skipGains = true; options.verbose = verbose;
[w] = ppr(f, g, h2q(h), 1/eta, degree, options);

end
