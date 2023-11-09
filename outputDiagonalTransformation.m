function [sigmaSquared, Tod] = outputDiagonalTransformation(v, w, Tin, Sigma, degree, verbose)
%outputDiagonalTransformation Compute the output-diagonal transformation for a polynomial control-affine dynamical system.
%
%   Usage: [sigmaSquared,Tbar] = outputDiagonalTransformation(v, w, Tin, Sigma, degree)
%
%   Inputs:
%       v,w     - cell arrays containing the polynomial energy function
%                 coefficients; these should already be in input-normal form.
%       Tin     - cell array containing the input-normal transformation
%                 coefficients.
%       Sigma   - diagonal matrix of the Hankel singular values of the
%                 linearized system.
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
%       Tod          - cell array containing the output-diagonal
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
if (nargin < 6)
    verbose = false;
    if (nargin < 5)
        degree = length(v) - 1;
    end
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

% if (degree < 1)
%     error('outputDiagonalTransformation: degree must be at least 1')
% end

% if (dv < degree + 1 || dw < degree + 1)
%     error('outputDiagonalTransformation: we need degree %d terms in the energy function', ...
%         degree + 1)
% end

% preallocate storage for the output T.
Tod = cell(1, degree); Tod{1} = speye(n);

[vtilde, wtilde] = transformEnergyFunctions(v, w, Tin, true); % Input-normal

for k = 3:degree+1
    if true
        [Nk] = equivalenceClassIndices(n, k);
        Nkhat = Nk(n + 1:end, :);
        
        M1 = Nk;
        M2Tr = zeros(n^k, size(Nkhat,1));
        for i=1:size(Nkhat,1)
            M2Tr(:,i) = vec(Sigma.^2  * reshape(Nkhat(i,:).',n,[]));
        end
        
        tempSize = size(Nk,1);
        Ainv = spdiags(1./sum(Nk,2),0,tempSize,tempSize); %diag(1./sum(Nk,2));
        SchurAInv = spdiags(1./diag(M2Tr.'*(speye(n^k) - Nk.'*Ainv*Nk) * M2Tr), 0, tempSize-n,tempSize-n);
        B = Nk * M2Tr; C = B.';
        
        M = 2*[M1; M2Tr.'];
        MMTinv = 1/4*[Ainv + Ainv*B*SchurAInv*C*Ainv , -Ainv*B*SchurAInv ;
            -SchurAInv*C*Ainv, SchurAInv];
        
        MPseudoInv = M.'*MMTinv;
        
        RHS = zeros(size(MPseudoInv,2),1);
        % Construct the RHS vector
        RHS = [];
        temp = zeros(size(Nk, 2), 1);
        for i = 2:k - 2
            j = k - i;
            temp = temp + vec(Tod{j}.' * Tod{i});
        end
        RHS = [RHS; -Nk * temp];
        temp = zeros(size(Nk(n + 1:end, :), 2), 1);
        for i = 2:k - 2
            j = k - i;
            temp = temp + vec(Tod{j}.' * Sigma .^ 2 * Tod{i});
        end
        for i = 3:k
            temp = temp + calTTv(Tod, i, k, wtilde{i}); % TODO: Accelerate this
        end
        RHS = [RHS; -Nk(n + 1:end, :) * temp];
        
        % TODO: Could cut down size by using kronecker repeated entries
                              
        %% Solve equations
        fprintf("Computing degree %i transformation... ", k-1)
        tic
        Tod{k - 1} = reshape(MPseudoInv * RHS,n,[]);             % Method 2
        fprintf("%f seconds. \n", toc)
    else
        [Nk] = equivalenceClassIndices(n, k);
        
        %$ Form input-normal equations coefficient matrix
        CoeffMatrix = 2 * Nk;
        
        % Construct the RHS vector
        RHS = [];
        temp = zeros(size(Nk, 2), 1);
        for i = 2:k - 2
            j = k - i;
            temp = temp + vec(Tod{j}.' * Tod{i});
        end
        RHS = [RHS; -Nk * temp];
        
        % \section
        %% Form output-diagonal equations coefficient matrix
        CoeffMatrix = [CoeffMatrix; 2 * Nk(n + 1:end, :) * kron(speye(n ^ (k - 1)), Sigma .^ 2)];
        
        temp = zeros(size(Nk(n + 1:end, :), 2), 1);
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
        
        %% Solve equations
        fprintf("Computing degree %i transformation... ", k-1)
        disp(size(indices))
        tic
        Tod{k - 1} = zeros(n, n ^ (k - 1));
        %     if k > 4
        %         Tod{k - 1}(indices) = CoeffMatrix \ RHS;                   % Method 1
        %     else
        
        Tod{k - 1}(indices) = pinv(full(CoeffMatrix)) * RHS;             % Method 2
        
        %         nullDir = rand(diff(size(CoeffMatrix)),1);                   % Method 3
        %         Tod{k - 1}(indices) = CoeffMatrix \ RHS  + null(CoeffMatrix) * nullDir;
        %     end
        
        % Symmetrize (needed when inverting the jacobian or something)
        %     for i=1:n
        %         Tod{k - 1}(i,:) = kronMonomialSymmetrize(Tod{k - 1}(i,:), n, k-1);
        %     end
        
        fprintf("%f seconds. \n", toc)
    end
end

%% Pluck out the singular value function coefficients

[vbar, wbar] = transformEnergyFunctions(vtilde, wtilde, Tod, true);

sigmaSquared = zeros(n, degree);

for k = 2:degree+1
    if verbose
        [N] = equivalenceClassIndices(n, k);
        
        fprintf("  - The largest entry in v%i is %.1e; ", k, max(abs(N * vbar{k}))) % Should be zero, other than the first time which is one
        fprintf("the largest off-diagonal entry in w%i is %.1e\n", k, max(abs(N(n + 1:end, :) * wbar{k}))) % Should be diagonal
        
        %         sigmaSquared{k - 1} = N(1:n, :) * wbar{k}; % Since the index set is already computed
        sigmaSquared(:,k - 1) = N(1:n, :) * wbar{k}; % Since the index set is already computed
    else
        indexSet = linspace(1, n ^ k, n);
        %         sigmaSquared{k - 1} = wbar{k}(indexSet);
        sigmaSquared(:,k - 1) = wbar{k}(indexSet);
    end
end

end

function [N] = equivalenceClassIndices(n, k)
%equivalenceClassIndices - For a k-way tensor of dimension n (n^k entries), compute the matrix N which combines the equivalence class entries. This matrix essentially goes from Kronecker product form to the unique monomial form.
%
% Usage: [N] = equivalenceClassIndices(n,k)
%
% TODO: currently not very efficient

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
% N = (repmat(linclassidx.', length(unique_values), 1) == unique_values); % Binary equivalence class matrix
N = sparse((repmat(linclassidx.', length(unique_values), 1) == unique_values)); % Binary equivalence class matrix

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
