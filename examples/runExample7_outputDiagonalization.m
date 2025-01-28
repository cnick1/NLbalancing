function [sigmaSquared] = runExample7_outputDiagonalization()
%runExample7_outputDiagonalization Runs 3D aircraft stall model output
% diagonalization.
%
%   Usage:  [sigmaSquared] = runExample7_outputDiagonalization()
%
%   runExample7_outputDiagonalization() runs the aircraft stall
%   stabilization model from Garrard 1977 [1]. This example specifically
%   tests the outputDiagonalTransformation() function.
%
%   Outputs:
%       sigmaSquared - coefficients of the square of the singular value
%                      functions
%
%   Reference: [1] W. L. Garrard and J. M. Jordan, “Design of nonlinear
%               automatic flight control systems,” Automatica, vol. 13,
%               no. 5, pp. 497–505, Sep. 1977,
%               doi: 10.1016/0005-1098(77)90070-x
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 7\n')
eta = 1; % values should be between -\infty and 1.

[f, g, h] = getSystem7();

%  Compute the energy functions
degree = 8;

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

%% Display the singular value functions

fprintf("\n  - Squared singular value functions:\n\n")

syms z
for i = 1:n
    fprintf("         𝜎_%i^2(z) = ", i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 3))
end

end
