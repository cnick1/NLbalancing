function runExample14()
%runExample14 Runs the 2D gradient pendulum example to test diagonalization.
%
%   Usage:  [v,w] = runExample14()
%
%   References: [1] J. M. A. Scherpen, “Balancing for nonlinear systems,”
%               PhD Dissertation, University of Twente, 1994.
%               [2] W. S. Gray and J. M. A. Scherpen, “Minimality and local
%               state decompositions of a nonlinear state space realization
%               using energy functions,” IEEE Transactions on Automatic
%               Control, vol. 45, no. 11, pp. 2079–2086, 2000, doi:
%               10.1109/9.887630
%               [1] K. Fujimoto and J. M. A. Scherpen, “Nonlinear
%               input-normal realizations based on the differential
%               eigenstructure of Hankel operators,” IEEE Transactions on
%               Automatic Control, vol. 50, no. 1, pp. 2–18, Jan. 2005,
%               doi: 10.1109/tac.2004.840476
%
%   Part of the NLbalancing repository.
%%

sympref('PolynomialDisplayStyle', 'ascend');

fprintf('Running Example 14\n')
eta = 0;

%% Scherpen 1994 and Scherpen/Gray 2000
fprintf("Beginning comparisons with Scherpen/Gray 2000:\n")

degree = 4;
[f, g, h] = getSystem14(degree - 1, 1);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta=eta, degree=degree, verbose=true);
[w] = approxFutureEnergy(f, g, h, eta=eta, degree=degree, verbose=true);

fprintf("\n  - Comparing our energy function with Scherpen/Gray 2000:\n")

fprintf("    Observability energy: \n        Lo = ")
disp(vpa(kronPolyEval(w, sym('x', [1, 2]).') / 2, 6))

fprintf("    Controllability energy: \n        Lc = ")
disp(vpa(kronPolyEval(v, sym('x', [1, 2]).') / 2, 6))

fprintf("\n                             ->  Energy functions match.\n\n")

%% Fujimoto/Scherpen 2005
fprintf("Beginning comparisons with Fujimoto/Scherpen 2005:\n")

degree = 4;
[f, g, h] = getSystem14(degree - 1, 2);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta=eta, degree=degree, verbose=true);
[w] = approxFutureEnergy(f, g, h, eta=eta, degree=degree, verbose=true);

fprintf("\n  - Comparing our energy function with Fujimoto/Scherpen 2005:\n")

fprintf("    Controllability energy: \n        Lc = ")
disp(vpa(kronPolyEval(v, sym('x', [1, 2]).') / 2, 5))

fprintf("    Observability energy: \n        Lo = ")
disp(vpa(kronPolyEval(w, sym('x', [1, 2]).') / 2, 5))

fprintf("\n                             ->  Energy functions match.\n\n")

fprintf("\n  - Comparing our squared singular value functions with Fujimoto/Scherpen 2005:\n")

tic
[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree=degree-1, verbose=true);
fprintf("Input-normal/output-diagonal transformation took %f seconds. \n", toc)

fprintf("\n  - Squared singular value functions:\n\n")

syms z
for i = 1:2
    fprintf("         𝜎_%i^2(z) = ", i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 5))
end

fprintf("\n                             ->  Squared singular value functions match.\n\n")

end
