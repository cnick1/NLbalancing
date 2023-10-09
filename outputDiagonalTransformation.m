function [singularValueFunSquared, Tod] = outputDiagonalTransformation(v, w, Tin, Sigma, degree, verbose)
%  Computes a polynomial approximation to the output-diagonal transformation
%  for a system with polynomial nonlinearities as
%
%     [singularValueFun,Tbar] = outputDiagonalTransformation(v,w,degree)
%
%  using polynomial approximations to the past and future energy functions.
%  The variables v and w contain the coefficients to the past and future energy
%  functions, respectively, already in input-normal form. Hence v2 should
%  be identity and v3 and on should be zero. Terms out to v{degree+1} and
%  w{degree+1} must be defined in the input.  Thus,
%
%      E_past(x) = 1/2 (kron(x,x))
%                = 0.5*kronPolyEval(v,x,degree+1)
%
%  and
%
%      E_future(x) = 1/2 ( w{2}kron(x,x) + ... + w{degree+1}kron(kron...,x),x) )
%                  = 0.5*kronPolyEval(w,x,degree+1)
%
%  The balancing transformation then has the form
%
%      x = T{1}z + T{2}kron(z,z) + ... + T{degree}kron(kron...,z),z)
%
%  where E_past(x) = 1/2 ( z.'z ) and E_future(x) = 1/2 ( z.'diag(sigma^2)z )
%  in the z coordinates.  The singular value functions are sigma.
%
%
%  Author: Nick Corbin, UCSD
%
%  License: MIT
%
%  Part of the NLbalancing repository.
%%
if (nargin < 4)
    verbose = false;
end

if (verbose)
    fprintf('Computing the degree %d output-diagonal balancing transformation\n', degree)
end

% Create a vec function for readability
vec = @(X) X(:);

validateattributes(v, {'cell'}, {})
validateattributes(w, {'cell'}, {})
dv = length(v);
dw = length(w);
n = sqrt(length(v{2}));

if (degree < 1)
    error('outputDiagonalTransformation: degree must be at least 1')
end

% if (dv < degree + 1 || dw < degree + 1)
%     error('outputDiagonalTransformation: we need degree %d terms in the energy function', ...
%         degree + 1)
% end

% preallocate storage for the output T.
Tod = cell(1, degree); Tod{1} = speye(n);

[vtilde, wtilde] = transformEnergyFunctions(v, w, Tin); % Input-normal

for k = 3:degree
    [Nk] = equivalenceClassIndices(n, k);

    % \section{}
    %$ Form input-normal equations coefficient matrix
    CoeffMatrix = 2 * Nk;

    % Construct the RHS vector
    RHS = [];
    temp = sparse(size(Nk, 2), 1);
    for i = 2:k - 2
        j = k - i;
        temp = temp + vec(Tod{j}.' * Tod{i});
    end
    RHS = [RHS; -Nk * temp];

    % \section
    %% Form output-diagonal equations coefficient matrix
    CoeffMatrix = [CoeffMatrix; 2 * Nk(n + 1:end, :) * kron(speye(n ^ (k - 1)), Sigma .^ 2)];

    temp = sparse(size(Nk(n + 1:end, :), 2), 1);
    for i = 2:k - 2
        j = k - i;
        temp = temp + vec(Tod{j}.' * Sigma .^ 2 * Tod{i});
    end
    for i = 3:k
        temp = temp + calTTv(Tod, i, k, wtilde{i}); % TODO: Accelerate this
    end
    RHS = [RHS; -Nk(n + 1:end, :) * temp];

    %% Form `flexibility' equations (Kronecker product repeated entries
    [linclassidx] = referenceElementMap(n, k - 1);

    linclassidx(linclassidx) = []; % Basically remove the reference element so one is nonzero and the rest we eliminate
    idxs = vec((n * (linclassidx - 1) + (1:n)).');

    %% Set extra parameter equation
    % TODO: find best parameters; for now just solve "a" solution with mldivide
    % parameterEqsRHS = 1;
    % parameterEqsCoeff = zeros(1,n^k);
    % parameterEqsCoeff(4) = 1;

    % CoeffMatrix(:,4) = [];
    % indices(4) = [];

    % idxs = [idxs; 16]; % alternatively could solve minimum norm solution

    %% Assemble equations
    CoeffMatrix(:, idxs) = []; % Kronecker flexibility

    % Form index set for the nonzero transformation components
    indices = 1:n ^ k; indices(idxs) = [];

    Tod{k - 1} = sparse(n, n ^ (k - 1));
    Tod{k - 1}(indices) = CoeffMatrix \ RHS;

end

[vbar, wbar] = transformEnergyFunctions(vtilde, wtilde, Tod);

singularValueFunSquared = cell(1, degree);

for k = 2:degree
    if verbose
        [N] = equivalenceClassIndices(n, k);

        fprintf("The largest entry in v%i is %.1e; ", k, max(abs(N * vbar{k}))) % Should be zero, other than the first time which is one
        fprintf("the largest off-diagonal entry in w%i is %.1e\n", k, max(abs(N(n + 1:end, :) * wbar{k}))) % Should be diagonal

        singularValueFunSquared{k - 1} = N(1:n, :) * wbar{k}; % Since the index set is already computed
    else
        indexSet = linspace(1, n ^ k, n);
        singularValueFunSquared{k - 1} = wbar{k}(indexSet);
    end
end

end

function [vtilde, wtilde] = transformEnergyFunctions(v, w, T)
%transformEnergyFunctions - Given a transformation T, compute the transformed energy functions given by the coefficients vtilde, wtilde.
%
% Syntax: [vtilde,wtilde] = transformEnergyFunctions(v,w,T)
%
% TODO: option to skip input if input normal transformation

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

function [N] = equivalenceClassIndices(n, k)
%equivalenceClassIndices - For a k-way tensor of dimension n (n^k entries), compute the matrix N which combines the equivalence class entries. This matrix essentially goes from Kronecker product form to the unique monomial form.
%
% Syntax: [N] = equivalenceClassIndices(n,k)
%

%% Compute equivalence class index sets
% We will compute a vector like [1 5 5 7 5 7 7 8]; the number in the vector
% corresponds to the reference element for the equivalence class, and all
% of the entries with the same number are in the same equivalence class.
% For an n-dimensional k-order tensor (n^k entries), there are
% nchoosek(n+k-1,k) unique entries (distinct equivalence classes, i.e.
% monomials)
[linclassidx] = referenceElementMap(n, k);

% We can use find(linclassidx == linclassidx(6)) to find all the indices
% that are equivalent to the index 6 for example
% find(linclassidx == linclassidx(2))

%% Form input-normal equations
% Get unique values from the input vector
[unique_values, ~, uidx] = unique(linclassidx, 'stable');

% Count the occurrences of each unique value, for sorting
multiplicity = accumarray(uidx, 1);
[~, sortMonomials] = sort(multiplicity); % Normalize to one permutation, i.e. reference element
unique_values = unique_values(sortMonomials);

% Construct the coefficient matrix
N = (repmat(linclassidx.', length(unique_values), 1) == unique_values); % Binary equivalence class matrix

end

function [linclassidx] = referenceElementMap(n, k)
%referenceElementMap - For a k-way tensor of dimension n (n^k entries),
%compute the vector that maps each element to its reference element.
%
% Syntax: [linclassidx] = referenceElementMap(n, k)
%

%% Compute equivalence class index sets
% We will compute a vector like [1 5 5 7 5 7 7 8]; the number in the vector
% corresponds to the reference element for the equivalence class, and all
% of the entries with the same number are in the same equivalence class.
% For an n-dimensional k-order tensor (n^k entries), there are
% nchoosek(n+k-1,k) unique entries (distinct equivalence classes, i.e.
% monomials)

% Construct matrix ind where each row is the multi-index for one element of X
idx = tt_ind2sub(ones(1, k) * n, (1:n ^ k)');

% Find reference index for every element in the tensor - this is to its
% index in the symmetrized tensor. This puts every element into a 'class'
% of entries that will be the same under symmetry.
classidx = sort(idx, 2); % Normalize to one permutation, i.e. reference element
mult = [1 cumprod(ones(1, k - 1) * n)]; % Form shifts
linclassidx = (classidx - 1) * mult' + 1; % Form vector that maps to the reference elements

end
