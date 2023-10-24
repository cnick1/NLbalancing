function runExample12(degree)
%runExample12 Runs the example to test diagonalization. The system is a
%   polynomial approximation of the 2D model from Fujimoto and Scherpen
%   2001, 2005, 2010 [1-3].
%
%   Usage:  [v,w] = runExample12(degree)
%
%   Inputs:
%       degree - desired degree of the energy function approximation
%
%   Outputs:
%
%   The value of eta is set below.
%
%   References: [1] K. Fujimoto and J. M. A. Scherpen, “Model reduction
%                for nonlinear systems based on the differential
%                eigenstructure of Hankel operators,” in Proceedings of
%                the 40th IEEE Conference on Decision and Control (Cat.
%                No.01CH37228), IEEE, 2001. doi: 10.1109/cdc.2001.980322
%               [2] K. Fujimoto and J. M. A. Scherpen, “Nonlinear
%                input-normal realizations based on the differential
%                eigenstructure of Hankel operators,” IEEE Transactions
%                on Automatic Control, vol. 50, no. 1, pp. 2–18, Jan.
%                2005, doi: 10.1109/tac.2004.840476
%               [3] K. Fujimoto and J. M. A. Scherpen, “Balanced
%                realization and model order reduction for nonlinear
%                systems based on singular value analysis,” SIAM Journal
%                on Control and Optimization, vol. 48, no. 7, pp.
%                4591–4623, Jan. 2010, doi: 10.1137/070695332
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 12\n')
eta = 0;

if nargin < 1
    degree = 8;
end

[f, g, h] = getSystem12(degree - 1);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, true);
[w] = approxFutureEnergy(f, g, h, eta, degree, true);

%% Compute the input-normal transformation approximation
[sigma, Tin] = inputNormalTransformation(v, w, degree - 1, false);

%% Compute the output-diagonal transformation approximation, also giving the squared singular value functions
[sigmaSquared, Tod] = outputDiagonalTransformation(v, w, Tin, diag(sigma), degree, true);

%% Plot the singular value functions
n = length(f{1});
z = linspace(- .5, .5, 51);
figure; hold on;
for i = 1:n
    plot(z, sqrt(polyval(flip(sigmaSquared(i, :)), z)))
end

end
