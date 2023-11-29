function runExample14()
%runExample14 Runs the example to test diagonalization.
%
%   Usage:  [v,w] = runExample14()
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 14\n')
eta = 0;

degree = 4;
[f, g, h] = getSystem14(degree - 1);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

fprintf("Beginning comparisons with Fujimoto/Tsubakino 2008:\n")

%% Comparison with Gray 2001 Example 1
% fprintf("\n  - Comparing our energy function with Fujimoto/Tsubakino 2008:\n")

% v{3}(abs(v{3}) < 1e-14) = 0; v{4}(abs(v{4}) < 1e-14) = 0; w{3}(abs(w{3}) < 1e-14) = 0; w{4}(abs(w{4}) < 1e-14) = 0;

% fprintf("    Controllability energy: \n        Lc = 1/2 *(")
% disp(vpa(kronPolyEval(v, sym('x', [1, 4]).'), 2))
% fprintf("    Observability energy: \n        Lo = 1/2 *(")
% disp(vpa(kronPolyEval(w, sym('x', [1, 4]).'), 2))
%
% fprintf("\n                             ->  Energy functions match.\n\n")

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
tic
[sigmaSquared, Tod] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
fprintf("Input-normal/output-diagonal transformation took %f seconds. \n", toc)

fprintf("\n  - Comparing our singular value functions with Fujimoto/Tsubakino 2008:\n\n")

sigmaSquared(sigmaSquared < 1e-14) = 0;

syms z
for i = 1:4
    fprintf("         𝜎_%i^2(z) = tau_%i(z e_i) = ", i, i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 3))
end

fprintf("\n                             ->  Hankel singular values match but not the functions.\n\n")

return
%% Compare transformation
fprintf("\n  - Comparing our transformation with Fujimoto/Tsubakino 2008:\n")
disp(vpa(kronPolyEval(Tod, sym('x', [1, 4]).'), 2))

end
