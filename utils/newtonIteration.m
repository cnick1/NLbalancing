function z = newtonIteration(x, f, J)
%newtonIteration Solve ƒ(z) = x for z using Newton iteration
%
%   Usage:  z = newtonIteration(x, f, J)
%
%   Inputs:     x - point at which to evaluate the transformation
%               f - function handle for the function ƒ(z)
%               J - function handle for the Jacobian J(z) = ∂ƒ(z)/∂z
%
%   Outputs:    z - the value of z = ƒ⁻¹(x)
%
%   Description: Assume we have a function x = ƒ(z); given z, we can
%   obtain x by evaluating the function forwards. However, what if we are
%   given x and wish to find the z such that ƒ(z) = x? This is the inverse
%   problem of finding the function z = ƒ⁻¹(x); in general, this is
%   nontrivial for nonlinear functions ƒ(z).
%
%   One method for solving this problem is the famous Newton iteration
%   [1]. It is straightforward to derive Newton's method, but for space I
%   will omit the derivation and simply give the formula. Define the
%   function g(z) = ƒ(z) - x, where x is fixed so not a variable; the goal
%   is to find roots of the equation g(z) = 0. Then, given an initial guess
%   z_0, we can iteratively improve our approximation as
%
%       z_(n+1)  =  z_n - g(z_n) / g'(z_n)
%
%   where g'(z_n) is the derivative of the function (assuming we can
%   evaluate it as well). The extension to multivariate polynomials is
%   straightforward: now z and g(z) are vectors, so g'(z) is the Jacobian
%   matrix (matrix of partial derivatives ∂g_i(z)/∂z_j), which since the
%   vector x is fixed is the same as the Jacobian of f(z). Defining the
%   Jacobian as J(z) = ∂g(z)/∂z, the Newton iteration formula then becomes
%
%       z_(n+1)  =  z_n - J(z_n) \ g(z_n)
%
%   which is the solution of a linear system of equations at each step.
%   This is done in line 55. Note that in the code, we do not store z_n,
%   we overwrite it at each step.
%
%   Specific usage: For usage with input-normal/output-diagonal
%   transformation x=Φ(z), the function handles would be
%       f = @(z) kronPolyEval(TinOd, z)
%       J = @(z) jcbn(TinOd, z)
%   This can all be done in one line if desired:
%       z0 = newtonIteration(x0, @(z) kronPolyEval(TinOd, z), @(z) jcbn(TinOd, z));
%
%   For usage with balancing transformation x=Φbar(zbar), we need to
%   compose the transformations as x=Φ(𝝋(zbar))
%       f = @(z) PhiBar(z,TinOd,sigmaSquared)
%       J = @(z) PhiBarJacobian(z,TinOd,sigmaSquared)
%   This can all be done in one line if desired:
%       z0 = newtonIteration(x0, @(z) PhiBar(z,TinOd,sigmaSquared), @(z) PhiBarJacobian(z,TinOd,sigmaSquared));
%
%   References: [1] https://en.wikipedia.org/wiki/Newton's_method
%                   (see the multidimensonal formulation section)
%
%   Part of the NLbalancing repository.
%%

% Solve for z0 initial condition with a Newton type iteration
fprintf("    Using Newton iteration to find transformed initial condition ... ")

maxIter = 10; % can change max iterations
z = x.*0; % Default initial guess of z=0
tol = 1e-14; % can edit tolerance

for iter = 1:maxIter
    % Update guess iteratively
    z = z - J(z)\(f(z)-x); % Proper Newton iteration
    
    % Break upon convergence within tolerance
    if norm(f(z)-x) < tol
        fprintf("converged in %i iterations. ",iter)
        break;
    end
end

fprintf(['\n         -> Initial condition: z0 = [', repmat('%2.2e ', 1, numel(z)), '], '], z);
fprintf('       error: %2.2e \n', norm(f(z)-x))

end