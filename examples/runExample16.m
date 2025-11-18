function runExample16()
%runExample16 Runs the 6D triple pendulum example to test diagonalization.
%
%   Usage:  runExample16()
%
%   References: [1] A. J. Krener, “Reduced order modeling of nonlinear
%               control systems,” in Analysis and Design of Nonlinear
%               Control Systems, Springer Berlin Heidelberg, 2008, pp.
%               41–62. doi: 10.1007/978-3-540-74358-3_4
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 16\n')
eta = 0;

degree = 4;
[f, g, h] = getSystem16(degree - 1);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta=eta, degree=degree);
[w] = approxFutureEnergy(f, g, h, eta=eta, degree=degree);;

fprintf("Beginning comparisons with Krener 2008:\n")

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
tic
[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree=degree-1, verbose=true);
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
disp(vpa(kronPolyEval(TinOd, sym('x', [1, 6]).'), 2))

end
