function [fbal,gbal,hbal,Tbal,sigmaSquared] = getBalancedRealization(f,g,h,nvp)
%getBalancedRealization Returns a balanced realization for f,g,h in Kronecker polynomial form
%
%   Usage: [fbal,gbal,hbal,Tbal,sigmaSquared] = getBalancedRealization(f,g,h)
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
%        degree - desired degree of the balanced realization. Ex. for linear
%                 dynamics, choosing degree=1 will produce quadratic
%                 approximations for the energy functions, linear approximations
%                 for the transformations, and the output will be a linear
%                 balanced realization. The default will be the degree of the
%                 drift term f(x) in the dynamics, i.e. we will balance
%                 everything available, no more and no less.
%        transformationDegree - desired degree of balancing transformation.
%                 Generally, this should be determined by the degree of
%                 balanced realization desired, i.e. the degree of the
%                 transformation through the option "degree". However, this
%                 additional option permits the user to override and for
%                 example use a linear transformation while producing a
%                 nonlinear transformed model. The most likely use case
%                 therefore corresponds to using transformationDegree=1 to
%                 compute quadratic energy functions and a linear
%                 transformation.
%       verbose - optional argument to print runtime information
%
%   Output:
%       fbal,gbal,hbal - cell arrays containing the polynomial coefficients for
%                        the drift, input, and output.
%       Tbal           - cell array containing the polynomial coefficients
%                        for the balancing transformation.
%       sigmaSquared   - squared singular value functions
%
%   Description: For a control affine dynamical system, we wish to compute a
%   balanced realization also in control-affine form via the balancing
%   transformation  x = ̅Φ(z̄):
%            ẋ = f(x) + g(x) u,      ->       ż̄ = f̃(z̄) + g̃(z̄) u
%            y = h(x).               ->       y = h̃(z̄)
%   The balanced realization is computed in three steps using helper functions:
%       1) Compute the balancing energy functions
%       2) Compute the balancing transformation
%       3) Compute the realization of the dynamics in the balanced coordinates
%   Additional details on each step can be found in the respective sections,
%   with further details in the helper functions.
%
%   References: [1] J. M. A. Scherpen, “Balancing for nonlinear systems,”
%               Systems & Control Letters, vol. 21, no. 2, pp. 143–153, Aug.
%               1993, doi: 10.1016/0167-6911(93)90117-o.
%               [2] J. M. A. Scherpen and A. J. Van Der Schaft, “Normalized
%               coprime factorizations and balancing for unstable nonlinear
%               systems,” International Journal of Control, vol. 60, no. 6, pp.
%               1193–1222, Aug. 1993, doi: 10.1080/00207179408921517.
%               [3] J. M. A. Scherpen, “Balancing for nonlinear systems,”
%               University of Twente, 1994.
%               [4] K. Fujimoto and J. M. A. Scherpen, “Balanced realization and
%               model order reduction for nonlinear systems based on singular
%               value analysis,” SIAM Journal on Control and Optimization, vol.
%               48, no. 7, pp. 4591–4623, Jan. 2010, doi: 10.1137/070695332.
%               [5] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” IEEE
%               Transactions on Automatic Control, vol. 70, no. 5, pp.
%               3088–3100, May 2025, doi: 10.1109/tac.2024.3494472.
%               [6] N. A. Corbin, A. Sarkar, J. M. A. Scherpen, and B. Kramer,
%               “Scalable computation of input-normal/output-diagonal balanced
%               realization for control-affine polynomial systems,” Systems &
%               Control Letters, vol. 204, p. 106178, Oct. 2025, doi:
%               10.1016/j.sysconle.2025.106178.
%
%   Part of the NLbalancing repository.
%
%   See also: approxPastEnergy, approxFutureEnergy, balancingTransformation, transformDynamics
arguments
    f
    g
    h
    nvp.degree = length(f)
    nvp.transformationDegree = []
    nvp.eta = 0
    nvp.verbose = false
end
if isempty(nvp.transformationDegree)
    nvp.transformationDegree = nvp.degree;
end
n = length(f{1});

if nvp.verbose
    fprintf("    Drift Dynamics:\n      ")
    dispKronPoly(f)
end

%% Step 1) Compute the energy functions
% The H∞ past and future energy functions (generalizations of controllability
% and observability) are defined as
%           𝓔⁻(x) = minᵤ J(x,u) = ½∫ η||y||² + ||u||²   dt
%           𝓔⁺(x) = minᵤ J(x,u) = ½∫  ||y||² + ||u||²/η dt
% where -1 < η < ∞.  Open-loop balancing is achieved with η = 0, and closed-loop
% balancing is achieved with η = 1. The functions approxPastEnergy() and
% approxFutureEnergy() compute polynomial approximations to the energy functions
% in Kronecker product form.
[v] = approxPastEnergy(f, g, h, eta=nvp.eta, degree=nvp.transformationDegree+1, verbose=nvp.verbose);
[w] = approxFutureEnergy(f, g, h, eta=nvp.eta, degree=nvp.transformationDegree+1, verbose=nvp.verbose);

if nvp.verbose
    fprintf("    Energy functions:\n      ")
    dispKronPoly(v,n=n), fprintf("\b")
    dispKronPoly(w,n=n)
end

%% Step 2) Compute the balancing transformation
% The nonlinear balancing transformation x = ̅Φ(z̄) = Φ(𝝋(z̄)) is the
% composition of the input-normal/output-diagonal transformation x = Φ(z) and
% the scaling transformation z = 𝝋(z̄). The balancing transformation puts the
% energy functions in the diagonal balanced form
%           𝓔⁻( ̅Φ(z̄)) = 1/2 z̄ᵀ Σ⁻¹(z̄) z̄
%           𝓔⁺( ̅Φ(z̄)) = 1/2 z̄ᵀ  Σ(z̄)  z̄
% where Σ(z̄) is the diagonal matrix of singular value functions ̅σᵢ(z̄ᵢ). The
% function balancingTransformation() computes an approximate polynomial
% expansion for the balancing transformation.
[Tbal, sigmaSquared] = balancingTransformation(v, w, degree=nvp.degree, verbose=nvp.verbose);

if nvp.verbose
    fprintf("    Balancing transformation:\n      ")
    dispKronPoly(Tbal)
end

%% Step 3) Compute the transformed dynamics
% Given the balancing transformation x = ̅Φ(z̄, now we seek to represent the
% dynamics for the control-affine system
%        ẋ = f(x) + g(x) u, y = h(x)
% in the new coordinates as
%        ż̄ = f̃(z̄) + g̃(z̄) u, y = h̃(z̄)
% The function transformDynamics() computes an approximate polynomial expansion
% for the transformed dynamics without inverting the nonlinear Jacobian.
[fbal,gbal,hbal] = transformDynamics(f,g,h,Tbal,degree=nvp.degree);

if nvp.verbose
    fprintf("    Transformed Drift Dynamics:\n      ")
    dispKronPoly(fbal)
end

end