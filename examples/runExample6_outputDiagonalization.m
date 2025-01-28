function [sigmaSquared] = runExample6_outputDiagonalization()
%runExample6_outputDiagonalization Runs the finite element beam example to test diagonalization.
%
%   Usage:  [sigmaSquared] = runExample6_outputDiagonalization()
%
%   Outputs:
%       sigmaSquared - coefficients of the square of the singular value
%                      functions
%
%   The value of eta is set below.
%
%   Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” 2023.
%
%   Part of the NLbalancing repository.
%%

numEls = 2;

eta = 1;

fprintf('Running Example 6\n')

[f, g, h] = getSystem6(numEls, 1);
g = g(1:2);
scaleFactor = 1e-9;
f = cellfun(@(x) x * scaleFactor, f, 'un', 0);

% A =  f{1}; B = g{1}; C = h{1};

% Compute the energy functions
degree = 6;

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
    plot(z, sqrt(polyval(flip(sigmaSquared(i, :)), z)))
end

fprintf("\n  - Squared singular value functions:\n\n")

syms z
for i = 1:n
    fprintf("         𝜎_%i^2(z) = ", i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 3))
end

end
