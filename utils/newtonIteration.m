function z = newtonIteration(x, g, J, maxIter, verbose)
%newtonIteration Solve g(z) = x for z using Newton iteration
%
%   Usage:  z = newtonIteration(x, g, J)
%
%   Inputs:     x - point at which to evaluate the transformation
%               g - function handle for the function g(z)
%               J - function handle for the Jacobian J(z) = ∂g(z)/∂z
%         maxIter - maximum number of iterations
%         verbose - boolean, whether to print runtime info
%
%   Outputs:    z - the value of z = g⁻¹(x)
%
%   Description: Assume we have a function x = g(z); given z, we can
%   obtain x by evaluating the function forwards. However, what if we are
%   given x and wish to find the z such that g(z) = x? This is the inverse
%   problem of finding the function z = g⁻¹(x); in general, this is
%   nontrivial for nonlinear functions g(z).
%
%   One method for solving this problem is the famous Newton iteration
%   [1]. It is straightforward to derive Newton's method, but for space I
%   will omit the derivation and simply give the formula. Define the
%   function ƒ(z) = g(z) - x, where x is fixed so not a variable; the goal
%   is to find roots of the equation ƒ(z) = 0. Then, given an initial guess
%   z_0, we can iteratively improve our approximation as
%
%       zₙ₊₁  =  zₙ - ƒ(zₙ) / ƒ'(zₙ)
%
%   where ƒ'(zₙ) is the derivative of the function (assuming we can
%   evaluate it as well). The extension to multivariate polynomials is
%   straightforward: now z and ƒ(z) are vectors, so ƒ'(z) is the Jacobian
%   matrix (matrix of partial derivatives ∂gᵢ(z)/∂zⱼ), which since the
%   vector x is fixed is the same as the Jacobian of ƒ(z). Defining the
%   Jacobian as J(z) = ∂ƒ(z)/∂z, the Newton iteration formula then becomes
%
%       zₙ₊₁  =  zₙ - J(zₙ)⁻¹ ƒ(zₙ)
%
%   which is the solution of a linear system of equations at each step.
%   This is done in line 79. Note that in the code, we do not store zₙ,
%   we overwrite it at each step.
%
%   Specific usage: For usage with input-normal/output-diagonal
%   transformation x = Φ(z), the function handles would be
%       g = @(z) kronPolyEval(TinOd, z)
%       J = @(z) jcbn(TinOd, z)
%   This can all be done in one line if desired:
%       z0 = newtonIteration(x0, @(z) kronPolyEval(TinOd, z), @(z) jcbn(TinOd, z));
%
%   For usage with balancing transformation x = ̅Φ(z̄), we need to
%   compose the transformations as x = Φ(𝝋(z̄))
%       g = @(z) PhiBar(z,TinOd,sigmaSquared)
%       J = @(z) PhiBarJacobian(z,TinOd,sigmaSquared)
%   This can all be done in one line if desired:
%       z0 = newtonIteration(x0, @(z) PhiBar(z,TinOd,sigmaSquared), @(z) PhiBarJacobian(z,TinOd,sigmaSquared));
%
%   References: [1] https://en.wikipedia.org/wiki/Newton's_method
%                   (see the multidimensonal formulation section)
%
%   Part of the NLbalancing repository.
%%
if nargin < 5
    verbose = false;
    if nargin < 4
        maxIter = 10;
        if nargin < 3
            error("Must supply x, f, and J")
        end
    end
end

% Solve for z0 initial condition with a Newton type iteration
if verbose; fprintf("    Using Newton iteration to find transformed initial condition ... "); end

z = x.*0; % Default initial guess of z=0
tol = 1e-14; % can edit tolerance

for iter = 1:maxIter
    % Update guess iteratively
    z = z - J(z)\(g(z)-x); % Proper Newton iteration
    
    % Break upon convergence within tolerance
    if norm(g(z)-x) < tol
        if verbose; fprintf("converged in %i iterations. ",iter); end
        break;
    end
end

if verbose; fprintf(['\n         -> Initial condition: z0 = [', repmat('%2.2e ', 1, numel(z)), '], '], z); end
if verbose; fprintf('       error: %2.2e \n', norm(g(z)-x)); end

end