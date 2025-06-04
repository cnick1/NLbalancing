function [sigmaSquared, TinOd, vbar, wbar] = inputNormalOutputDiagonalTransformation(v, w, degree, verbose)
%inputNormalOutputDiagonalTransformation Compute the input-normal/output-diagonal transformation x = Φ(z) for a polynomial control-affine dynamical system.
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
%                      NOT given by sigmaSquared.^(½). Use the function
%                      utils/polySqrt() to get the coefficients of the
%                      singular value functions.
%
%       TinOd        - cell array containing the output-diagonal
%                      transformation coefficients.
%
%   Description: We compute a transformation x = Φ(z) that makes the energy
%   functions input-normal
%           𝓔⁻(Φ(z)) = ½ zᵀz
%   and output-diagonal
%           𝓔⁺(Φ(z)) = ½ zᵀ Σ²(z) z,
%   where Σ²(z) is the diagonal matrix of squared singular value functions
%   σᵢ²(zᵢ). In our case, the energy functions are computed as polynomials
%   in the original coordinates:
%           𝓔⁻(x) = ½ ( v₂ᵀ(z⊗z) + v₃ᵀ(z⊗z⊗z) + ... ),
%           𝓔⁺(x) = ½ ( w₂ᵀ(z⊗z) + w₃ᵀ(z⊗z⊗z) + ... ).
%   In general, in transformed coordinates, the coefficients will be
%           𝓔⁻(Φ(z)) = ½ ( ṽ₂ᵀ(z⊗z) + ṽ₃ᵀ(z⊗z⊗z) + ... )
%           𝓔⁺(Φ(z)) = ½ ( w̃₂ᵀ(z⊗z) + w̃₃ᵀ(z⊗z⊗z) + ... )
%   Input-normal corresponds to ṽ₂ being identity and ṽ₃ and above being
%   zero. Output-diagonal corresponds to w̃₂ being a diagonal matrix and
%   w̃₃ and above being diagonal tensors. The transformation is computed by
%   representing it as
%            x = Φ(z)
%              = T₁z + T₂(z⊗z) + ... + Td(z...⊗z)
%   and deriving the conditions on T₁, T₂, etc. to ensure the desired
%   structure in the transformed coefficients. Terms out to v{degree+1} and
%   w{degree+1} must be defined in the input to the function.
%
%   The approach taken here is a two-step approach which enables sparsity.
%   First, the linear input-normal/output-diagonal transformation is
%   applied. Then, the nonlinear transformation components are computed in
%   the transformed coordinates. The two transformations are then combined.
%   Additional details can be found in [1].
%
%   References: [1] N. A. Corbin, A. Sarkar, J. M. A. Scherpen, and B. Kramer,
%                “Scalable computation of input-normal/output-diagonal balanced
%                realization for control-affine polynomial systems,” Oct. 2024,
%                doi: 10.48550/arXiv.2410.22435
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
    fprintf('Computing the degree %d input-normal/output-diagonal balancing transformation...\n', degree)
end

% Create a vec function for readability
vec = @(X) X(:);

validateattributes(v, {'cell'}, {})
validateattributes(w, {'cell'}, {})
dv = length(v); %#ok<NASGU>
dw = length(w); %#ok<NASGU>
n = sqrt(numel(v{2}));

method = 'two-step'; % 'two-step' is best
switch method
    case 'BorisJeffSerkan'
        %% Boris/Jeff/Serkan's method, as described in the Part 2 paper
        % In the part 2 paper, the transformation they compute corresponds
        % to the out-dated "smooth" singular value function definition
        % rather than the "axis" singular value function definition.
        % Practically speaking, this means they only truly satisfy
        % input-normal and not output-diagonal. In the details of how the
        % method works, this amounts to computing the squared singular value
        % functions by picking off the diagonals and throw away rest. In
        % theory this is fine to get the squared singular value functions,
        % but the transformation that gets you there is only input-normal
        % and NOT output-diagonal (hence the need to throw away the
        % off-diagonal terms, as they are nonzero)
        %% Compute the input-normal transformation approximation
        [~, Tin] = inputNormalTransformation(v, w, degree, [], false);
        
        TinOd = Tin;
        
    case 'direct'
        %% Our first approach: input normal, then output diagonal
        % Here, we compute first the input-normal transformation from the
        % Part 2 paper above, but then we compute an output diagonal
        % transformation those transformed coordinates. This is a
        % deprecated approach, as we figured out that we can do only the
        % linear term in the first transformation and the nonlinear
        % transformation terms only need to be computed once.
        %% Compute the input-normal transformation approximation
        [sigma, Tin] = inputNormalTransformation(v, w, degree, [], false);
        
        %% Compute the output-diagonal transformation approximation, also giving the squared singular value functions
        [~, Tod] = outputDiagonalTransformation(v, w, Tin, diag(sigma), degree, false);
        
        %% Combine the transformations
        % This is approximate; some information lost due to truncation of
        % the product of the two polynomial transformations. Two-step
        % approach below is better.
        TinOd = composeTransformations(Tin, Tod);
        
    case 'two-step'
        %% Two-step Input-Normal/Output-Diagonal Transformation
        % This is the approach described in Corollary 1 of [1]. First we
        % compute the linear transformation, then we compute the nonlinear
        % terms in the transformed coordinates, which enables sparsity, and
        % finally we combine the linear and nonlinear transformations.
        
        %% Step 1: Compute the linear input-normal/output-diagonal transformation
        % In principle, v{2} and w{2} would be vectors; I have rewritten
        % the code to use the factoredMatrix class for V2 and W2 though, so
        % v{2} and w{2} are actually factoredMatrix objects and reshape
        % does nothing.
        V2 = reshape(v{2}, n, n); W2 = reshape(w{2}, n, n);
        
        if isa(V2, 'factoredMatrixInverse') && isa(W2, 'factoredMatrix')
            % handle factored matrix class
            Rinv = cholinv(V2);  % V₂ = "P⁻¹" = (R⁻¹*R⁻ᵀ)⁻¹ = R*Rᵀ
            L = chol(W2);        % W₂ = "Q" = L*Lᵀ
            % [~, Xi, V] = svd(L.' * Rinv.'); % from Theorem 2
            [V, Xi, ~] = svd(Rinv * L); % just to avoid transposing things
            
            % Truncate transformation to degree r
            %   Xi = Xi(1:r,1:r); V = V(:,1:r);
            
            Tin = Rinv.' * V; % linear input-normal/output-diagonal transformation
        else
            % this should never be called now by my implementation; can
            % probably get rid of it
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
            % Truncate transformation to degree r
            %   Xi = Xi(1:r,1:r); V = V(:,1:r);
            
            Tin = R.' \ V; % linear input-normal/output-diagonal transformation
        end
        
        % Transform the energy functions using Tin
        [vtilde, wtilde] = transformEnergyFunctionsLinear(v, w, Tin);
        
        % Name Ṽ₂ and W̃₂; in principle they would be the first two
        % entries, but analytically we know what they are
        % V2tilde = reshape(vtilde{2},n,n); W2tilde = reshape(wtilde{2},n,n);
        V2tilde = speye(n); vtilde{2} = vec(V2tilde);
        W2tilde = sparse(diag(diag(Xi))) .^ 2; wtilde{2} = vec(W2tilde);
        
        %% Step 2: Compute the higher-order terms in the second transformation
        % Preallocate the cell array, the first term is identity
        Tod = cell(1, degree);
        Tod{1} = speye(n);
        
        % Compute the higher-order terms according to Corollary 1 [1]
        for k = 3:degree + 1
            
            [Nk, Nkhat] = equivalenceClassIndices(n, k);
            
            %% Form input-normal equations coefficient matrix
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
            fprintf("    Computing degree %i coefficient... ", k - 1)
            %         disp(size(indices))
            tic
            Tod{k - 1} = zeros(n, n ^ (k - 1));
            %     if k > 4
            %             spparms('spumoni',2)
            Tod{k - 1}(indices) = CoeffMatrix \ RHS;                 % Method 1: matlab uses sparse QR from SuiteSparseQR
            % Tod{k - 1}(indices) = lsqminnorm(CoeffMatrix, RHS);      % Method 2: minimum norm solution
            
            %% Optional: adjust transformation by a random element from the null-space
            %             nullDir = rand(diff(size(CoeffMatrix)),1);
            %             Tod{k - 1}(indices) = Tod{k - 1}(indices).' + null(full(CoeffMatrix)) * nullDir;
            
            fprintf("completed in %f seconds. \n", toc)
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
        
        fprintf("      - The largest entry in v%i is %.1e; ", k, max(abs(N * vbar{k}))) % Should be zero, other than the first time which is one
        fprintf("the largest off-diagonal entry in w%i is %.1e\n", k, max(abs(N(n + 1:end, :) * wbar{k}))) % Should be diagonal
        
        %         sigmaSquared{k - 1} = N(1:n, :) * wbar{k}; % Since the index set is already computed
        sigmaSquared(:, k - 1) = N(1:n, :) * wbar{k}; % Since the index set is already computed
    else
        indexSet = linspace(1, n ^ k, n);
        %         sigmaSquared{k - 1} = wbar{k}(indexSet);
        sigmaSquared(:, k - 1) = wbar{k}(indexSet);
    end
end

if verbose
    % Plot the squared singular value functions
    z = linspace(- 1, 1, 101);
    figure; hold on; title("Singular value functions")
    for i = 1:n
        plot(z, real(sqrt(polyval(flip(sigmaSquared(i, :)), z))))
        % plot(z, polyval(flip(sigmaSquared(i, :)), z)) % plot squared singular value functions
    end
    set(gca,'yscale','log')
    xlabel('z_i','Interpreter','TeX'); ylabel('\sigma_i','Interpreter','TeX'); legend('\sigma_1','\sigma_2','Interpreter','TeX')
end

end

function [N, Nhat] = equivalenceClassIndices(n, k)
%equivalenceClassIndices Compute the equivalence class entry mapping
% For a k-way tensor of dimension n (which has n^k entries), compute the
% matrix N which combines the equivalence class entries. This matrix
% essentially maps from Kronecker product form to the unique monomial form.

%% Compute equivalence class index sets
[linclassidx] = referenceElementMap(n, k);

%% Form input-normal equations
% Get unique values from the input vector
diagIdxs = linspace(1, n^k, n); offdiagIdxs = setdiff(1:n^k, diagIdxs);
linclassidx(diagIdxs) = [];
[~, ~, uidx] = unique(linclassidx, 'stable');
Nhat = sparse(uidx, offdiagIdxs, 1, nchoosek(n+k-1, k) - n, n^k);
N = [sparse(1:n, linspace(1, n^k, n), 1); Nhat];

end

function [linclassidx] = referenceElementMap(n, k)
%referenceElementMap Compute the reference element mapping for an equivalence class of monomials
% For a k-way tensor of dimension n (which has n^k entries), compute the
% vector that maps each element to its reference element.

%% Compute equivalence class index sets
% We will compute a vector like [1 5 5 7 5 7 7 8]; the number in the vector
% corresponds to the reference element for the equivalence class, and all
% of the entries with the same number are in the same equivalence class.
% For an n-dimensional k-order tensor (n^k entries), there are nchoosek(n+k-1,k)
% unique entries (distinct equivalence classes, i.e. monomials)

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
%       vtilde,wtilde - cell arrays containing the transformed polynomial
%                       energy function coefficients
%
%   Description: Given a linear transformation x=Tz, compute the transformed energy functions
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
