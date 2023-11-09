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
%   Reference: [0] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” 2023.
%              [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki, “Nonlinear
%               balanced truncation: Part 1—computing energy functions,” arXiv,
%               Dec. 2022. doi: 10.48550/ARXIV.2209.07645
%              [2] B. Kramer, S. Gugercin, and J. Borggaard, “Nonlinear balanced
%               truncation: Part 2—model reduction on manifolds,” arXiv, Feb. 2023.
%               doi: 10.48550/ARXIV.2302.02036
%              [3] J. Borggaard and L. Zietsman, “On approximating polynomial-
%               -quadratic regulator problems,” IFAC-PapersOnLine, vol. 54, no. 9,
%               pp. 329–334, 2021, doi: 10.1016/j.ifacol.2021.06.090
%
%  Part of the NLbalancing repository.
%%

n = 8;
eta = 0.9;

[f, g, h, zInit] = getSystem3(n, 4, 4, 0.001, 0);

fprintf('Running Example 3\n')

% Compute the energy functions
degree = 5;

[v] = approxPastEnergy(f, g, h, eta, degree, true);
[w] = approxFutureEnergy(f, g, h, eta, degree, true);

%% Compute the input-normal transformation approximation
[sigma, Tin] = inputNormalTransformation(v, w, degree - 1, false);

%% Compute the output-diagonal transformation approximation, also giving the squared singular value functions
[sigmaSquared, Tod] = outputDiagonalTransformation(v, w, Tin, diag(sigma), degree - 1, true);

%% Plot the singular value functions
n = length(f{1});
z = linspace(- .5, .5, 51);
figure; hold on;
for i = 1:n
%     plot(z, sqrt(polyval(flip(sigmaSquared(i, :)), z)))
    plot(z, (polyval(flip(sigmaSquared(i, :)), z)))
end

fprintf("\n  - Singular value functions:\n\n")

syms z
for i=1:n
    fprintf("         𝜎_%i^2(z) = ",i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 3))
end

end
