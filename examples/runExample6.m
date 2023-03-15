function [v, w] = runExample6(degree)
%runExample6 Runs the example from the paper
%
%   Usage:  [v,w] = runExample6(degree)
%
%   where
%         degree          is the degree of energy function approximations
%
%         v,w             are coefficients of the past and future energy
%                         function approximations, respectively.
%
%   The value of eta is set below.
%
%   Part of the NLbalancing repository.
%%

fprintf('Running Example 6\n')

eta = 0.1; % values should be between -\infty and 1.
% eta=0.1 corresponds to gamma= 1.0541...
% since eta = 1 - 1/gamma^2;

fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

if (nargin < 1)
  degree = 4;
end

[A, B, C, N, Q, f, g, h] = getSystem6();

%  Compute the polynomial approximations to the energy functions
[v] = approxPastEnergy(A, N, g(1:2), C, eta, degree, true);
[w] = approxFutureEnergy(A, N, g(1:2), C, eta, degree, true);

end
