function [sigmaSquared] = runExample3_outputDiagonalization()
%runExample3_outputDiagonalization Runs the Burgers discretization to test diagonalization.
%
%   Usage:  [sigmaSquared] = runExample6_outputDiagonalization()
%
%   Outputs:
%       sigmaSquared - coefficients of the square of the singular value
%                      functions
%
%   The value of eta is set below.
%
%   Reference: [0] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗∞
%               energy functions for polynomial control-affine systems,"
%               IEEE Transactions on Automatic Control, pp. 1–13, 2024,
%               doi: 10.1109/tac.2024.3494472
%              [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki,
%               “Scalable computation of energy functions for nonlinear
%               balanced truncation,” Computer Methods in Applied Mechanics
%               and Engineering, vol. 427, p. 117011, Jul. 2024, doi:
%               10.1016/j.cma.2024.117011
%              [2] B. Kramer, S. Gugercin, and J. Borggaard, “Nonlinear balanced
%               truncation: Part 2—model reduction on manifolds,” arXiv, Feb. 2023.
%               doi: 10.48550/arXiv.2302.02036
%              [3] J. Borggaard and L. Zietsman, “On approximating polynomial-
%               -quadratic regulator problems,” IFAC-PapersOnLine, vol. 54, no. 9,
%               pp. 329–334, 2021, doi: 10.1016/j.ifacol.2021.06.090
%
%  Part of the NLbalancing repository.
%%

n = 16;
eta = 0.9;

[f, g, h, zInit] = getSystem3(n, 4, 4, 0.001, 0);

fprintf('Running Example 3\n')

% Compute the energy functions
degree = 4;

[v] = approxPastEnergy(f, g, h, eta, degree, true);
[w] = approxFutureEnergy(f, g, h, eta, degree, true);

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
tic
[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
fprintf("Input-normal/output-diagonal transformation took %f seconds. \n", toc)

%% Plot the squared singular value functions
n = length(f{1});
z = linspace(- .5, .5, 51);
figure; hold on;
for i = 1:n
    %     plot(z, sqrt(polyval(flip(sigmaSquared(i, :)), z)))
    plot(z, (polyval(flip(sigmaSquared(i, :)), z)))
end

fprintf("\n  - Squared singular value functions:\n\n")

syms z
for i = 1:n
    fprintf("         𝜎_%i^2(z) = ", i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 3))
end

end
