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

method = 1; % Methods 1,2,4 work, 1 is fastest
for k = 3:degree + 1
    switch method
        case 5 % try trivial output diagonal solution
            [Nk] = equivalenceClassIndices(n, k);
            
            %% Form output-diagonal equations coefficient matrix
            CoeffMatrix = 2 * Nk(n + 1:end, :) * kron(speye(n ^ (k - 1)), Sigma .^ 2);
            
            temp = zeros(size(Nk(n + 1:end, :), 2), 1);
            for i = 2:k - 2
                j = k - i;
                temp = temp + vec(Tod{j}.' * Sigma .^ 2 * Tod{i});
            end
            for i = 3:k
                temp = temp + calTTv(Tod, i, k, wtilde{i}); % TODO: Accelerate this
            end
            RHS = -Nk(n + 1:end, :) * temp;
            
            %% Form `flexibility' equations (Kronecker product repeated entries
            [linclassidx] = referenceElementMap(n, k - 1);
            
            linclassidx(linclassidx) = []; % Basically remove the reference element so one is nonzero and the rest we eliminate
            idxs = vec((n * (linclassidx - 1) + (1:n)).');
            
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
            Tod{k - 1}(indices) = CoeffMatrix \ RHS;                   % Method 1
            
            fprintf("%f seconds. \n", toc)
        case 4 % Another try at pseudoinverse
            fprintf("Computing degree %i transformation... ", k - 1)
            tic
            Nk = equivalenceClassIndices(n, k);
            Nkhat = Nk(n + 1:end, :);
            
            % Form Minv
            M1 = full(2* Nk);
            
            % M2 Method 1: Kronecker identity
            tempMT = zeros(n ^ k, size(Nkhat, 1));
            for i = 1:size(Nkhat, 1)
                tempMT(:, i) = vec(Sigma .^ 2 * reshape(Nkhat(i, :).', n, []));
            end
            M2 = 2*(tempMT).'; clear tempMT
            
            MPseudoInv = [pinv(M1 - M1*(pinv(M2)*M2)), pinv(M2 - M2*(pinv(M1)*M1))];
            
            % Form RHS
            % Construct the RHS vector
            b1 = zeros(n^k, 1);
            b2 = zeros(n^k, 1);
            for i = 2:k - 2
                j = k - i;
                b1 = b1 + vec(Tod{j}.' * Tod{i});
                b2 = b2 + vec(Tod{j}.' * Sigma .^ 2 * Tod{i});
            end
            for i = 3:k
                b2 = b2 + calTTv(Tod, i, k, wtilde{i}); % TODO: Accelerate this
            end
            b1 = -Nk * b1; 
            b2 = -Nkhat * b2;
            
%             %%
%             idx = tt_ind2sub(ones(1, k) * n, (1:n ^ k)');
%             classidx = sort(idx, 2); % Normalize to one permutation, i.e. reference element
%             mult = [1 cumprod(ones(1, k - 1) * n)]; % Form shifts
%             linclassidx = (classidx - 1) * mult' + 1; % Form vector that maps to the reference elements
%             uniquelinclassidx = unique(linclassidx);
%             b11 = -accumarray(linclassidx,b1); b1 = b11(uniquelinclassidx);
% %             linclassidx(linspace(1,n^k,n)) = [];
%             b22 = -accumarray(linclassidx,b2); b2 = b22(uniquelinclassidx);  
%             b2(find(sum(uniquelinclassidx == linspace(1,n^k,n),2))) = [];

            %% Solve equations
            Tod{k - 1} = zeros(n, n ^ (k - 1));
            Tod{k - 1}(:) = MPseudoInv * [b1 ; b2]; 
            fprintf("%f seconds. \n", toc)
        case 3
            %% Preliminary tools
            idx = tt_ind2sub(ones(1, k) * n, (1:n ^ k)');
            classidx = sort(idx, 2); % Normalize to one permutation, i.e. reference element
            mult = [1 cumprod(ones(1, k - 1) * n)]; % Form shifts
            linclassidx = (classidx - 1) * mult' + 1; % Form vector that maps to the reference elements
            
            multiplicity = accumarray(linclassidx, 1);
            %             c = accumarray(linclassidx,c);
            
            monomialSize = nchoosek(n+k-1,k);
            
            fprintf("Computing degree %i transformation... ", k - 1)
            tic
            
            %%
            b1 = zeros(n^k, 1);
            for i = 2:k - 2
                j = k - i;
                b1 = b1 + vec(Tod{j}.' * Tod{i});
            end
            NNb1 = kronMonomialSymmetrize(b1,n,k);
            
            
            % Method 2: Kronecker identity (new, still needs Nkhat, which needs Nk)
            Nk = equivalenceClassIndices(n, k);
            Nkhat = Nk(n + 1:end, :);
            M2Tr = zeros(n ^ k, size(Nkhat, 1));
            for i = 1:size(Nkhat, 1)
                M2Tr(:, i) = vec(Sigma .^ 2 * reshape(Nkhat(i, :).', n, []));
            end
            M2Tr = 2*sparse(M2Tr); M2 = M2Tr.';
            SchurAInv = M2 * (M2Tr - kronMonomialSymmetrize(M2,3,3).');
            
            % Method 1: Kronecker identity (old, needs Nk and Nkhat)
            %             M2Tr = zeros(n ^ k, size(Nkhat, 1));
            %             for i = 1:size(Nkhat, 1)
            %                 M2Tr(:, i) = vec(Sigma .^ 2 * reshape(Nkhat(i, :).', n, []));
            %             end
            %             M2Tr = sparse(M2Tr);
            %             Ainv = spdiags(1 ./ nonzeros(multiplicity), 0, monomialSize, monomialSize); %diag(1./sum(Nk,2));
            %             SchurAInv = spdiags(1 ./ diag(M2Tr.' * (speye(n ^ k) - Nk.' * Ainv * Nk) * M2Tr), 0, tempSize - n, tempSize - n);
            
            
            %%
            
            B = Nk * M2Tr; C = B.';
            
            M = 2 * [M1; M2];
            MMTinv = 1/4 * [Ainv + Ainv * B * SchurAInv * C * Ainv, -Ainv * B * SchurAInv;
                -SchurAInv * C * Ainv, SchurAInv];
            
            MPseudoInv = M.' * MMTinv;
            
            RHS = zeros(size(MPseudoInv, 2), 1);
            % Construct the RHS vector
            RHS = [b1; ];
            b2 = zeros(size(Nk(n + 1:end, :), 2), 1);
            for i = 2:k - 2
                j = k - i;
                b2 = b2 + vec(Tod{j}.' * Sigma .^ 2 * Tod{i});
            end
            for i = 3:k
                b2 = b2 + calTTv(Tod, i, k, wtilde{i}); % TODO: Accelerate this
            end
            RHS = [b1; -Nk(n + 1:end, :) * b2];
            
            % TODO: Could cut down size by using kronecker repeated entries
            
            %% Solve equations
            Tod{k - 1} = reshape(MPseudoInv * RHS, n, []); % Method 2
            fprintf("%f seconds. \n", toc)
        case 2 % First try at pseudoinverse construction
            fprintf("Computing degree %i transformation... ", k - 1)
            tic
            Nk = equivalenceClassIndices(n, k);
            Nkhat = Nk(n + 1:end, :);
            %             [Nk, Nkhat] = equivalenceClassIndices(n, k);
            
            M1 = Nk;
            % Method 1: Kronecker identity
%             M2Tr = sparse(n ^ k, size(Nkhat, 1));
%             for i = 1:size(Nkhat, 1)
%                 M2Tr(:, i) = vec(Sigma .^ 2 * reshape(Nkhat(i, :).', n, []));
%             end
            M2Tr = 2 * kron(speye(n ^ (k - 1)), Sigma .^ 2) * Nkhat.';
            % Method 2: reshaping and columnwise multiplication
            %         Nkhat_3D = reshape(permute(reshape(full(Nkhat), n, [], size(Nkhat, 1)), [1 3 2]),[],n);
            
            tempSize = size(Nk, 1);
            Ainv = spdiags(1 ./ sum(Nk, 2), 0, tempSize, tempSize); %diag(1./sum(Nk,2));
            SchurAInv = spdiags(1 ./ diag(M2Tr.' * (speye(n ^ k) - Nk.' * Ainv * Nk) * M2Tr), 0, tempSize - n, tempSize - n);
            B = Nk * M2Tr; C = B.';
            
            M = 2 * [M1; M2Tr.'];
            MMTinv = 1/4 * [Ainv + Ainv * B * SchurAInv * C * Ainv, -Ainv * B * SchurAInv;
                -SchurAInv * C * Ainv, SchurAInv];
            
            MPseudoInv = M.' * MMTinv;
            
            %             RHS = zeros(size(MPseudoInv, 2), 1);
            % Construct the RHS vector
            RHS = [];
            temp = zeros(size(Nk, 2), 1);
            for i = 2:k - 2
                j = k - i;
                temp = temp + vec(Tod{j}.' * Tod{i});
            end
            RHS = [RHS; -Nk * temp];
            temp = zeros(size(Nkhat, 2), 1);
            for i = 2:k - 2
                j = k - i;
                temp = temp + vec(Tod{j}.' * Sigma .^ 2 * Tod{i});
            end
            for i = 3:k
                temp = temp + calTTv(Tod, i, k, wtilde{i}); % TODO: Accelerate this
            end
            RHS = [RHS; -Nkhat * temp];
            
            % TODO: Could cut down size by using kronecker repeated entries
            
            %% Solve equations
            Tod{k - 1} = reshape(MPseudoInv * RHS, n, []); % Method 2
            fprintf("%f seconds. \n", toc)
            
        case 1
            [Nk, Nkhat] = equivalenceClassIndices(n, k);
%             Nkhat = Nk(n + 1:end, :);
%             Nkhat = Nk; Nkhat(linspace(1,n^k,n),:) = [];
            
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
            CoeffMatrix = [CoeffMatrix; 2 * Nkhat * kron(speye(n ^ (k - 1)), Sigma .^ 2)]; % TODO: kronecker rules
            
            temp = zeros(size(Nkhat, 2), 1);
            for i = 2:k - 2
                j = k - i;
                temp = temp + vec(Tod{j}.' * Sigma .^ 2 * Tod{i});
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
            Tod{k - 1}(indices) = CoeffMatrix \ RHS;                   % Method 1
            %     else
            
            %         Tod{k - 1}(indices) = pinv(full(CoeffMatrix)) * RHS; % Method 2
%             Tod{k - 1}(indices) = lsqminnorm(CoeffMatrix, RHS);
            
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

for k = 2:degree
    Tod{k} = sparse(kronMonomialUnsymmetrize(Tod{k}, n, k));
end

%% Pluck out the singular value function coefficients

[vbar, wbar] = transformEnergyFunctions(vtilde, wtilde, Tod, true);

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

% function [N] = ReducedNk(n, k)
% %ReducedNk - For a k-way tensor of dimension n (n^k entries), compute the matrix N which combines the equivalence class entries. This matrix essentially goes from Kronecker product form to the unique monomial form.
% %
% % Usage: [N] = ReducedNk(n,k)
% %
%
% referenceIndices = unique(referenceElementMap(n, k));
% m = length(referenceIndices);
%
% N = sparse(m, n ^ k); % initialize a matrix of size length(v) x n with zeros
%
% % Set the corresponding entries to 1
% N(1:m, referenceIndices) = speye(m);
%
% end
%
% function [N] = ReducedNkHat(n, k)
% %ReducedNkHat - For a k-way tensor of dimension n (n^k entries), compute the matrix N which combines the equivalence class entries. This matrix essentially goes from Kronecker product form to the unique monomial form.
% %
% % Usage: [N] = ReducedNkHat(n,k)
% %
%
% referenceIndices = referenceElementMap(n, k);
% diagonalIndices = linspace(1,n^k,n);
% referenceIndices(diagonalIndices) = [];
% referenceIndices = unique(referenceIndices);
%
% m = length(referenceIndices);
%
% N = sparse(m, n ^ k); % initialize a matrix of size length(v) x n with zeros
%
% % Set the corresponding entries to 1
% N(1:m, referenceIndices) = speye(m);
%
% end

function [N, Nhat] = equivalenceClassIndices(n, k)
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
% multiplicity = accumarray(linclassidx, 1);

% We can use find(linclassidx == linclassidx(6)) to find all the indices
% that are equivalent to the index 6 for example
% find(linclassidx == linclassidx(2))

%% Form input-normal equations
% Get unique values from the input vector
[~, ~, uidx] = unique(linclassidx, 'stable');

% % Method 2
% N = sparse(uidx,1:n^k,1);

diagIdxs = linspace(1,n^k,n); offdiagIdxs = setdiff(1:n^k, diagIdxs);
linclassidx(diagIdxs) = [];
[~, ~, uidx] = unique(linclassidx, 'stable');
Nhat = sparse(uidx,offdiagIdxs,1,nchoosek(n+k-1,k)-n,n^k);
N = [sparse(1:n, linspace(1,n^k,n), 1); Nhat];


% return
% Method 1
% Count the occurrences of each unique value, for sorting
% [unique_values, ~, uidx] = unique(linclassidx, 'stable');

% multiplicity = accumarray(uidx, 1);
% [~, sortMonomials] = sort(multiplicity); % Normalize to one permutation, i.e. reference element
% unique_values = unique_values(sortMonomials);
% N = sparse((repmat(linclassidx.', length(unique_values), 1) == unique_values)); % Binary equivalence class matrix
% Nhat = N(n+1:end,:);

% N = sparse((repmat(linclassidx.', length(unique_values), 1) == unique_values)); % Binary equivalence class matrix
% diagIdxs = linspace(1,n^k,n); 
% Nhat = N; Nhat(uidx(diagIdxs),:) = []; 

end

function [N, Nhat] = equivalenceClassIndicesNew(n, k)
%equivalenceClassIndices - For a k-way tensor of dimension n (n^k entries), compute the matrix N which combines the equivalence class entries. This matrix essentially goes from Kronecker product form to the unique monomial form.
%
% Usage: [N] = equivalenceClassIndices(n,k)
%
% TODO: currently not very efficient

%% Compute equivalence class index sets
[linclassidx] = referenceElementMap(n, k);

%% Method 2
% Get unique values from the input vector

temp(1:n^k, linclassidx) = speye(n^k);
N = temp(unique(linclassidx), linclassidx);

linclassidx2 = linclassidx;
linclassidx2(linspace(1,n^k,n)) = [];
Nhat = temp(unique(linclassidx2), linclassidx);


%% Method 1
% Get unique values from the input vector
% unique_values = unique(linclassidx);

% N = sparse((repmat(linclassidx.', length(unique_values), 1) == unique_values)); % Binary equivalence class matrix

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
