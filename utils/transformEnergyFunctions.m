function [vtilde, wtilde] = transformEnergyFunctions(v, w, T)
%transformEnergyFunctions Transforms the energy coefficients v and w by T.
%
%   Usage: [vtilde, wtilde] = transformEnergyFunctions(v, w, T)
%
%   Inputs:
%       v,w   - cell arrays containing the polynomial energy function coefficients
%       T     - cell array containing the polynomial transformation coefficients
%       TODO: option to skip input if input normal transformation
%
%   Output:
%       v,w   - cell arrays containing the transformed polynomial energy function
%               coefficients
%
%   Background: Given a transformation T, compute the transformed energy functions
%   given by the coefficients vtilde, wtilde. TODO: Add more details here.
%
%   Authors: Nick Corbin, UCSD
%

vec = @(X) X(:);

degree = length(w);
n = sqrt(length(v{2}));
V2 = reshape(v{2}, n, n);
W2 = reshape(w{2}, n, n);

vtilde = cell(1, degree);
wtilde = cell(1, degree);

vtilde{2} = vec(T{1}.' * V2 * T{1});
wtilde{2} = vec(T{1}.' * W2 * T{1});

for k = 3:degree
    vtilde{k} = vec(T{k - 1}.' * V2 * T{1}) + vec(T{1}.' * V2 * T{k - 1});
    wtilde{k} = vec(T{k - 1}.' * W2 * T{1}) + vec(T{1}.' * W2 * T{k - 1});

    for i = 2:k - 2
        j = k - i;
        vtilde{k} = vtilde{k} + vec(T{j}.' * V2 * T{i});
        wtilde{k} = wtilde{k} + vec(T{j}.' * W2 * T{i});
    end

    for i = 3:k
        vtilde{k} = vtilde{k} + calTTv(T, i, k, v{i});
        wtilde{k} = wtilde{k} + calTTv(T, i, k, w{i});
    end
end

end