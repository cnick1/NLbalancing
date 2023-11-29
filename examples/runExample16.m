function runExample16()
%runExample16 Runs the example to test diagonalization.
%
%   Usage:  [v,w] = runExample16()
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 16\n')
eta = 0;

degree = 4;
[f, g, h] = getSystem16(degree - 1);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

fprintf("Beginning comparisons with Krener 2008:\n")

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
tic
[sigmaSquared, Tod] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
fprintf("Input-normal/output-diagonal transformation took %f seconds. \n", toc)

fprintf("\n  - Comparing our Hankel singular values with Krener 2008:\n\n    ")

disp((sigmaSquared(:, 1).') .^ (1/2))
fprintf("                         ->  Hankel singular values match.\n\n")

fprintf("\n  - Comparing our squared singular value functions with Krener 2008:\n\n")

z = linspace(0, 1, 51);
figure; set(0, 'DefaultLineLineWidth', 2);
for i = 1:6
    semilogy(z, polyval(flip(sigmaSquared(i, :)), z)); hold on;
end
ylim([1e-5 1e3])
% semilogy(z, 2*polyval([2.5818e-5, 0, 2.2545e-07, 0, 8.2416e-05], z))

syms z
for i = 1:6
    fprintf("         𝜎_%i^2(z) = tau_%i(z e_i) = ", i, i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 3))
end

return
%% Compare transformation
fprintf("\n  - Comparing our transformation with Krener 2008:\n")
disp(vpa(kronPolyEval(Tod, sym('x', [1, 6]).'), 2))

end
