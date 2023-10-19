function runExample12()
%runExample12 Runs the example to test diagonalization. The system is a
%   polynomial approximation of the 2D model from Fujimoto and Scherpen
%   2001, 2005, 2010 [1-3].
%
%   Usage:  [v,w] = runExample12()
%
%   Inputs:
%       exportData      - Boolean variable to determine if
%                         plots/data are exported
%
%   Outputs:
%       v,w             - Coefficients of the past and future energy
%                         function approximations, respectively
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
fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

%  Compute the polynomial approximations to the future energy function
degree = 8;

[f, g, h] = getSystem12(degree)

[v] = approxPastEnergy(f, g, h, eta, degree + 1);
[w] = approxFutureEnergy(f, g, h, eta, degree + 1);

end
