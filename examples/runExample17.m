function sigmaSquared = runExample17(n)
%runExample17 Runs the example to test diagonalization.
%
%   Usage:  [v,w] = runExample17()
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 17\n')
eta = 0;

if nargin < 1
    n = 64;
end

degree = 4;
[f, g, h] = getSystem17(degree - 1, n / 2);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, true);
[w] = approxFutureEnergy(f, g, h, eta, degree, true);


%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
tic
[sigmaSquared, Tod] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
fprintf("Input-normal/output-diagonal transformation took %f seconds. \n", toc)

end
