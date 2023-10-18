function runExample12()
%runExample12 Runs the example to test diagonalization.
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
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 12\n')

eta = 0;
fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

%  Compute the polynomial approximations to the future energy function
degree = 8;

[f, g, h] = getSystem12(degree)

[v] = approxPastEnergy(f, f{2}, g, h, eta, degree+1);
[w] = approxFutureEnergy(f, f{2}, g, h, eta, degree+1);

end
