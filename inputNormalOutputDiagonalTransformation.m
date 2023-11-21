function [sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree, verbose)
%inputNormalOutputDiagonalTransformation Compute the input-normal/output-diagonal transformation for a polynomial control-affine dynamical system.
%
%   Usage: [sigmaSquared,Tbar] = outputDiagonalTransformation(v, w, Tin, Sigma, degree)
%
%   Inputs:
%       v,w     - cell arrays containing the polynomial energy function
%                 coefficients; these should already be in input-normal form.
%       degree  - desired degree of the computed transformation (default =
%                 degree of energy functions - 1).
%       verbose - optional argument to print runtime information.
%
%   Outputs:
%       sigmaSquared - an n×degree-1 matrix containing the coefficients of
%                      the square of the singular value functions. The
%                      first column corresponds to the square of the Hankel
%                      singular values, the next column corresponds to the
%                      degree 1 coefficients, etc. These can be plotted
%                      using polyval() (and flip()).
%
%                      Warning: Note that the singular value functions are
%                      NOT given by sigmaSquared.^(1/2). Use the function
%                      utils/polySqrt() to get the coefficients of the
%                      singular value functions.
%
%       TinOd        - cell array containing the output-diagonal
%                      transformation coefficients.
%
%   Background: Since the inputs should be in input normal form, v{2} should
%   be identity and v{3} and on should be zero. Terms out to v{degree+1} and
%   w{degree+1} must be defined in the input.  Thus,
%
%      E_past(x) = 1/2 (kron(x,x))
%                = 0.5*kronPolyEval(v,x,degree+1)
%      and E_future(x) = 1/2 ( w{2}kron(x,x) + ... + w{degree+1}kron(kron...,x),x) )
%                      = 0.5*kronPolyEval(w,x,degree+1)
%
%  The output-diagonal transformation then has the form
%
%      x = T{1}z + T{2}kron(z,z) + ... + T{degree}kron(kron...,z),z)
%
%  where E_past(x) = 1/2 ( z.'z ) and E_future(x) = 1/2 ( z.'diag(sigma^2)z )
%  in the z coordinates.  The singular value functions are sigma.
%
%  Author: Nick Corbin, UCSD
%
%  License: MIT
%
%  Part of the NLbalancing repository.
%%
if (nargin < 4)
    verbose = false;
    if (nargin < 3)
        degree = length(v) - 1;
    end
end

if (verbose)
    fprintf('Computing the degree %d input-normal/output-diagonal balancing transformation\n', degree)
end

% Create a vec function for readability
vec = @(X) X(:);

validateattributes(v, {'cell'}, {})
validateattributes(w, {'cell'}, {})
dv = length(v);
dw = length(w);
n = sqrt(length(v{2}));

method = 3; % 3 is best
switch method
    case 1 % Boris' method: input-normal only, pick off diagonals and throw away rest
        %% Compute the input-normal transformation approximation
        [~, Tin] = inputNormalTransformation(v, w, degree, [], false);
        
        TinOd = Tin;
        
    case 2 % Our first approach: input normal, then output diagonal
        %% Compute the input-normal transformation approximation
        [sigma, Tin] = inputNormalTransformation(v, w, degree, [], false);
        
        %% Compute the output-diagonal transformation approximation, also giving the squared singular value functions
        [~, Tod] = outputDiagonalTransformation(v, w, Tin, diag(sigma), degree, false);
        
        %% Combine the transformations
        % TODO: this is currently approximate; some information lost due to truncation as well
        TinOd = composeTransformations(Tin, Tod);
        
    case 3 % Our second approach: balance linear portion first, then balance nonlinear portions
        
        V2 = reshape(v{2}, n, n);
        W2 = reshape(w{2}, n, n);
        
        try
            R = chol(V2, 'lower'); % V2 = R*R.'
        catch
            warning("inputNormalOutputDiagonalTransformation: Cholesky factorizatin failed; trying sqrtm()")
            R = sqrtm(V2);
        end
        try
            L = chol(W2, 'lower'); % W2 = L*L.'
        catch
            warning("inputNormalOutputDiagonalTransformation: Cholesky factorizatin failed; trying sqrtm()")
            L = sqrtm(W2);
        end
        [~, Xi, V] = svd(L.' / R.'); % from Theorem 2
        
        %% Truncate transformation to degree r
        %   Xi = Xi(1:r,1:r);
        %   V = V(:,1:r);
        
        %% preallocate storage for the output T.
        
        Tin = R.' \ V;
        [vtilde, wtilde] = transformEnergyFunctionsLinear(v, w, Tin);
        V2tilde = speye(n); vtilde{2} = vec(V2tilde);
        % V2tilde = reshape(vtilde{2},n,n);
        W2tilde = sparse(diag(diag(Xi))) .^ 2; wtilde{2} = vec(W2tilde);
        % W2tilde = reshape(wtilde{2},n,n);
        
        Tod = cell(1, degree);
        
        Tod{1} = speye(n);
        
        for k = 3:degree + 1
            
            [Nk, Nkhat] = equivalenceClassIndices(n, k);
            
            %$ Form input-normal equations coefficient matrix
            CoeffMatrix = 2 * Nk;
            
            % Construct the RHS vector
            RHS = [];
            temp = zeros(size(Nk, 2), 1);
            for i = 2:k - 2
                j = k - i;
                temp = temp + vec(Tod{j}.' * V2tilde * Tod{i});
            end
            for i = 3:k
                temp = temp + calTTv(Tod, i, k, vtilde{i}); % TODO: Accelerate this
            end
            RHS = [RHS; -Nk * temp];
            
            %% Form output-diagonal equations coefficient matrix
            CoeffMatrix = [CoeffMatrix; 2 * Nkhat * kron(speye(n ^ (k - 1)), W2tilde)]; % TODO: kronecker rules
            
            temp = zeros(size(Nkhat, 2), 1);
            for i = 2:k - 2
                j = k - i;
                temp = temp + vec(Tod{j}.' * W2tilde * Tod{i});
            end
            for i = 3:k
                temp = temp + calTTv(Tod, i, k, wtilde{i}); % TODO: Accelerate this
            end
            RHS = [RHS; -Nkhat * temp];
            
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
            
            %% Solve equations
            fprintf("Computing degree %i transformation... ", k - 1)
            %         disp(size(indices))
            tic
            Tod{k - 1} = zeros(n, n ^ (k - 1));
            %     if k > 4
            Tod{k - 1}(indices) = CoeffMatrix \ RHS; % Method 1
            %     Tod{k - 1}(indices) = lsqminnorm(CoeffMatrix, RHS);      % Method 2
            
            fprintf("%f seconds. \n", toc)
        end
        
        %% Combine transformation with linear input normal transformation
        
        TinOd{1} = Tin;
        for k = 2:degree
            TinOd{k} = Tin * Tod{k};
        end
        
end

%% Symmetrize the transformation rows
for k = 2:degree
    %     Tod{k} = sparse(kronMonomialUnsymmetrize(Tod{k}, n, k));
    TinOd{k} = kronMonomialSymmetrize(TinOd{k}, n, k);
end

%% Pluck out the singular value function coefficients
[vbar, wbar] = transformEnergyFunctions(v, w, TinOd, true); % Could transform just the observability; could probably even just compute the diagonal entries

sigmaSquared = zeros(n, degree);

for k = 2:degree + 1
    if verbose
        [N] = equivalenceClassIndices(n, k);
        
        fprintf("  - The largest entry in v%i is %.1e; ", k, max(abs(N * vbar{k}))) % Should be zero, other than the first time which is one
        fprintf("the largest off-diagonal entry in w%i is %.1e\n", k, max(abs(N(n + 1:end, :) * wbar{k}))) % Should be diagonal
        
        %         sigmaSquared{k - 1} = N(1:n, :) * wbar{k}; % Since the index set is already computed
        sigmaSquared(:, k - 1) = N(1:n, :) * wbar{k}; % Since the index set is already computed
    else
        indexSet = linspace(1, n ^ k, n);
        %         sigmaSquared{k - 1} = wbar{k}(indexSet);
        sigmaSquared(:, k - 1) = wbar{k}(indexSet);
    end
end

end

function [N, Nhat] = equivalenceClassIndices(n, k)
%equivalenceClassIndices - For a k-way tensor of dimension n (n^k entries), compute the matrix N which combines the equivalence class entries. This matrix essentially goes from Kronecker product form to the unique monomial form.
%
% Usage: [N] = equivalenceClassIndices(n,k)
%

%% Compute equivalence class index sets
[linclassidx] = referenceElementMap(n, k);

%% Form input-normal equations
% Get unique values from the input vector

diagIdxs = linspace(1, n ^ k, n); offdiagIdxs = setdiff(1:n ^ k, diagIdxs);
linclassidx(diagIdxs) = [];
[~, ~, uidx] = unique(linclassidx, 'stable');
Nhat = sparse(uidx, offdiagIdxs, 1, nchoosek(n + k - 1, k) - n, n ^ k);
N = [sparse(1:n, linspace(1, n ^ k, n), 1); Nhat];

end

function [linclassidx] = referenceElementMap(n, k)
%referenceElementMap - For a k-way tensor of dimension n (n^k entries),
%compute the vector that maps each element to its reference element.
%
% Usage: [linclassidx] = referenceElementMap(n, k)
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

function [vtilde, wtilde] = transformEnergyFunctionsLinear(v, w, T)
%transformEnergyFunctionsLinear Transforms the energy coefficients v and w by T.
%
%   Usage: [vtilde, wtilde] = transformEnergyFunctions(v, w, T)
%
%   Inputs:
%       v,w         - cell arrays containing the polynomial energy function
%                     coefficients
%       T           - linear transformation coefficient
%
%   Output:
%       v,w   - cell arrays containing the transformed polynomial energy function
%               coefficients
%
%   Background: Given a linear transformation x=Tz, compute the transformed energy functions
%   given by the coefficients vtilde, wtilde. TODO: Add more details here.
%
%   Authors: Nick Corbin, UCSD
%

vec = @(X) X(:);

degree = length(w);
% n = sqrt(length(v{2}));
[n, r] = size(T);
V2 = reshape(v{2}, n, n);
W2 = reshape(w{2}, n, n);

vtilde = cell(1, degree);
wtilde = cell(1, degree);

vtilde{2} = vec(T.' * V2 * T);
wtilde{2} = vec(T.' * W2 * T);

for k = 3:degree
    vtilde{k} = calTTv({T}, k, k, v{k});
    wtilde{k} = calTTv({T}, k, k, w{k});
    
    vtilde{k} = kronMonomialSymmetrize(vtilde{k}, r, k);
    wtilde{k} = kronMonomialSymmetrize(wtilde{k}, r, k);
end

end
