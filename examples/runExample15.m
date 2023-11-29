function runExample15()
%runExample15 Runs the 4D double pendulum example to test diagonalization.
%
%   Usage:  [v,w] = runExample15()
%
%   References: [1]
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
[sigmaSquared, Tod] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
fprintf("Input-normal/output-diagonal transformation took %f seconds. \n", toc)

fprintf("\n  - Comparing our singular value functions with Fujimoto/Tsubakino 2008:\n\n")

% sigmaSquared(sigmaSquared < 1e-14) = 0;

syms z
for i = 1:4
    fprintf("         𝜎_%i^2(z) = tau_%i(z e_i) = ", i, i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 3))
end

fprintf("\n                             ->  Hankel singular values match but not the functions.\n\n")

end
