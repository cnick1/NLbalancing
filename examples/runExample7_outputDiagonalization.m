function [v, w] = runExample7_outputDiagonalization()
%runExample7_outputDiagonalization Runs 3D aircraft stall model output 
% diagonalization.
% 
%   Usage:  [v, w] = runExample7_outputDiagonalization()
%
%   runExample7_outputDiagonalization() runs the aircraft stall 
%   stabilization model from Garrard 1977 [1]. This example specifically
%   tests the outputDiagonalTransformation() function. 
%
%   Outputs:
%       v,w - coefficients of the past and future energy function
%             approximations, respectively.
%
%   Reference: [1] W. L. Garrard and J. M. Jordan, “Design of nonlinear
%               automatic flight control systems,” Automatica, vol. 13,
%               no. 5, pp. 497–505, Sep. 1977,
%               doi: 10.1016/0005-1098(77)90070-x
%
%   Part of the NLbalancing repository.
%%
vec = @(X) X(:);

[f, g, h] = getSystem7();
n = length(f{1});

fprintf('Running Example 7\n')
eta = 1; % values should be between -\infty and 1.

%  Compute the polynomial approximations to the future energy function
degree = 8;

[v] = approxPastEnergy(f, g, h, eta, degree, true);
[w] = approxFutureEnergy(f, g, h, eta, degree, true);

%% compute the input-normal transformation approximation
[sigma, Tin] = inputNormalTransformation(v, w, degree - 1, false);

%% compute the output-diagonal transformation approximation, also giving the squared singular value functions
[singularValueFunSquared, Tod] = outputDiagonalTransformation(v, w, Tin, diag(sigma), degree, true);

[c] = approximateSingularValueFunctions(Tin, w, sigma, degree - 2);

cell2mat(singularValueFunSquared)
[sigma, cell2mat(c)]

%%
[vtilde, wtilde] = transformEnergyFunctions(v, w, Tin); % Input-normal

wtilde{3}(1)



end
