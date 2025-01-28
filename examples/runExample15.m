function runExample15()
%runExample15 Runs the 4D double pendulum example to test diagonalization.
%
%   Usage:  [v,w] = runExample15()
%
%   References: [1] K. Fujimoto and D. Tsubakino, “On computation of
%               nonlinear balanced realization and model reduction,” in
%               2006 American Control Conference, IEEE, 2006. doi:
%               10.1109/acc.2006.1655399
%               [2] K. Fujimoto and D. Tsubakino, “Computation of nonlinear
%               balanced realization and model reduction based on Taylor
%               series expansion,” Systems & Control Letters, vol. 57, no.
%               4, pp. 283–289, Apr. 2008, doi: 10.1016/j.sysconle.2007.08.015
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 15\n')
eta = 0;

degree = 4;
[f, g, h] = getSystem15(degree - 1);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

fprintf("Beginning comparisons with Fujimoto/Tsubakino 2008:\n")

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
tic
[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
fprintf("Input-normal/output-diagonal transformation took %f seconds. \n", toc)

fprintf("\n  - Comparing our singular value functions with Fujimoto/Tsubakino 2008:\n\n")

% sigmaSquared(sigmaSquared < 1e-14) = 0;

syms z
for i = 1:4
    fprintf("         𝜎_%i^2(z) = ", i, i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 3))
end

fprintf("\n                             ->  Hankel singular values match but not the functions.\n\n")

end
