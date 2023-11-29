function runExample14_gradient()
%runExample14_gradient Runs the example to test diagonalization.
%
%   Usage:  [v,w] = runExample14_gradient()
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%

sympref('PolynomialDisplayStyle','ascend');

fprintf('Running Example 14\n')
eta = 0;

%% Scherpen 1994 and Scherpen/Gray 2000
fprintf("Beginning comparisons with Scherpen/Gray 2000:\n")

degree = 4;
[f, g, h] = getSystem14(degree - 1, 2);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, true);
[w] = approxFutureEnergy(f, g, h, eta, degree, true);

fprintf("\n  - Comparing our energy function with Scherpen/Gray 2000:\n")

fprintf("    Observability energy: \n        Lo = ")
disp(vpa(kronPolyEval(w, sym('x', [1, 2]).')/2, 6))

fprintf("    Controllability energy: \n        Lc = ")
disp(vpa(kronPolyEval(v, sym('x', [1, 2]).')/2, 6))

fprintf("\n                             ->  Energy functions match.\n\n")

%% Fujimoto/Scherpen 2005

fprintf("Beginning comparisons with Fujimoto/Scherpen 2005:\n")

degree = 4;
[f, g, h] = getSystem14(degree - 1, 3);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, true);
[w] = approxFutureEnergy(f, g, h, eta, degree, true);

fprintf("\n  - Comparing our energy function with Fujimoto/Scherpen 2005:\n")

fprintf("    Controllability energy: \n        Lc = ")
disp(vpa(kronPolyEval(v, sym('x', [1, 2]).')/2, 5))

fprintf("    Observability energy: \n        Lo = ")
disp(vpa(kronPolyEval(w, sym('x', [1, 2]).')/2, 5))

fprintf("\n                             ->  Energy functions match.\n\n")

fprintf("\n  - Comparing our squared singular value functions with Fujimoto/Scherpen 2005:\n")

tic
[sigmaSquared, Tod] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
fprintf("Input-normal/output-diagonal transformation took %f seconds. \n", toc)

fprintf("\n  - Squared singular value functions:\n\n")

syms z
for i = 1:2
    fprintf("         𝜎_%i^2(z) = ", i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 5))
end

fprintf("\n                             ->  Squared singular value functions match.\n\n")


end
