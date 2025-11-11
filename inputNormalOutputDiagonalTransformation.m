function [sigmaSquared, TinOd, vbar, wbar] = inputNormalOutputDiagonalTransformation(v, w, degree, verbose)
%inputNormalOutputDiagonalTransformation Return a polynomial input-normal/output-diagonal transformation x = Φ(z).
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
%                      NOT given by sigmaSquared.^(½).
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
%                realization for control-affine polynomial systems,” Systems &
%                Control Letters, vol. 204, p. 106178, Oct. 2025, doi:
%                10.1016/j.sysconle.2025.106178.

%
%  Author: Nick Corbin, UCSD
%
%  License: MIT
%
%  Part of the NLbalancing repository.
%%
arguments
    v cell
    w cell
    degree = length(v) - 1
    verbose = false
end
vec = @(X) X(:); % Create a vec function for readability
if verbose
    fprintf('Computing the degree %d input-normal/output-diagonal balancing transformation...\n', degree)
end

n = sqrt(numel(v{2}));

%% Two-step Input-Normal/Output-Diagonal Transformation
% This is the approach described in Corollary 1 of [1]. First we compute the
% linear transformation, then we compute the nonlinear terms in the transformed
% coordinates, which enables sparsity, and finally we combine the linear and
% nonlinear transformations.

%% Step 1: Compute the linear input-normal/output-diagonal transformation
% The code has been written to avoid computing v{2} and w{2} explicitly;
% they are stored as factoredMatrix objects, where we compute directly
% their Cholesky factors for square-root balancing. Here, we just retreive
% the square-root factors.
Rinv = cholinv(v{2});       % V₂ = "P⁻¹" = (R⁻ᵀ*R⁻¹)⁻¹ = R*Rᵀ
L = chol(w{2});             % W₂ = "Q"                 = L*Lᵀ
[V, Xi, U] = svd(Rinv * L); % same as UΣV=LᵀR⁻ᵀ from Theorem 2, just avoids transposing

% Truncate transformation to order r
% r = 2;
% Xi = Xi(1:r,1:r); U = U(:,1:r); V = V(:,1:r);

% Construct linear input-normal/output-diagonal transformation and inverse
Tin = invertibleMatrix(Rinv.'*V,  Xi\U.'*L.');

% Transform the energy functions using Tin
[vtilde, wtilde] = transformEnergyFunctionsLinear(v, w, Tin);

% Name Ṽ₂ and W̃₂; in principle they would be the first two entries, but
% analytically we know what they are
% V2tilde = reshape(vtilde{2},n,n); W2tilde = reshape(wtilde{2},n,n);
V2tilde = speye(n); vtilde{2} = vec(V2tilde);
W2tilde = sparse(diag(diag(Xi))) .^ 2; wtilde{2} = vec(W2tilde);

%% Step 2: Compute the higher-order terms in the second transformation
% Preallocate the cell array, the first term is identity
Tod = cell(1, degree);
Tod{1} = speye(n);

% Compute the higher-order terms according to Corollary 1 [1]
for k = 3:degree + 1
    fprintf("    Computing degree %i coefficient... ", k - 1); tic
    
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
    
    %% Form `flexibility' equations (Kronecker product repeated entries)
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
    indices = 1:n^k; indices(idxs) = [];
    
    %% Solve equations
    Tod{k-1} = zeros(n, n^(k-1));
    Tod{k-1}(indices) = CoeffMatrix \ RHS;                     % Method 1: matlab uses sparse QR from SuiteSparseQR
    % Tod{k - 1}(indices) = lsqminnorm(CoeffMatrix, RHS);      % Method 2: minimum norm solution
    
    fprintf("completed in %f seconds. \n", toc)
end

%% Combine transformation with linear input-normal transformation
TinOd = cell(1,degree);
TinOd{1} = Tin;
for k = 2:degree
    TinOd{k} = Tin * Tod{k};
    TinOd{k} = kronMonomialSymmetrize(TinOd{k}, n, k); % Symmetrize the transformation rows
end

%% Pluck out the singular value function coefficients
[vbar, wbar] = transformEnergyFunctions(v, w, TinOd, true); % Could transform just the observability; could probably even just compute the diagonal entries

sigmaSquared = zeros(n, degree);
for k = 2:degree + 1
    if verbose
        [N] = equivalenceClassIndices(n, k);
        
        fprintf("      - The largest entry in v%i is %.1e; ", k, max(abs(N * vbar{k}))) % Should be zero, other than the first time which is one
        fprintf("the largest off-diagonal entry in w%i is %.1e\n", k, max(abs(N(n + 1:end, :) * wbar{k}))) % Should be diagonal
        
        sigmaSquared(:, k - 1) = N(1:n, :) * wbar{k}; % Since the index set is already computed
    else
        indexSet = linspace(1, n ^ k, n);
        sigmaSquared(:, k - 1) = wbar{k}(indexSet);
    end
end

if verbose
    % Plot the squared singular value functions
    z = linspace(- 1, 1, 101);
    figure; hold on; title("Singular value functions")
    for i = 1:n
        plot(z, real(sqrt(polyval(flip(sigmaSquared(i, :)), z))))
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
classidx = sort(idx, 2);                  % Normalize to one permutation, i.e. reference element
mult = [1 cumprod(ones(1, k - 1) * n)];   % Form shifts
linclassidx = (classidx - 1) * mult' + 1; % Form vector that maps to the reference elements

end

function [vtilde, wtilde] = transformEnergyFunctionsLinear(v, w, T)
%transformEnergyFunctionsLinear Transforms the energy coefficients v and w by T.
%
%   Usage: [vtilde, wtilde] = transformEnergyFunctionsLinear(v, w, T)
%
%   Inputs:
%       v,w - cell arrays containing the polynomial energy function coefficients
%       T   - linear transformation coefficient
%
%   Output:
%       vtilde,wtilde - cell arrays containing the transformed polynomial
%                       energy function coefficients
%
%   Description: Consider past and future energy functions given by the
%   polynomial expansions
%           𝓔⁻(x) = ½ ( v₂ᵀ(z⊗z) + v₃ᵀ(z⊗z⊗z) + ... ),
%           𝓔⁺(x) = ½ ( w₂ᵀ(z⊗z) + w₃ᵀ(z⊗z⊗z) + ... ),
%   and consider a polynomial transformation
%            x = Φ(z)
%              = T₁z + T₂(z⊗z) + ... + Td(z...⊗z).
%   As shown in Lemma 1 in [1], the energy functions can be expressed in the
%   transformed z coordinates as
%           𝓔⁻(Φ(z)) = ½ ( ṽ₂ᵀ(z⊗z) + ṽ₃ᵀ(z⊗z⊗z) + ... )
%           𝓔⁺(Φ(z)) = ½ ( w̃₂ᵀ(z⊗z) + w̃₃ᵀ(z⊗z⊗z) + ... )
%   where the transformed coordinates are computed using the calligraphic T
%   notation according to
%                 ₖ                     ₖ
%           ṽₖᵀ = ∑ vⱼᵀ 𝓣ⱼ,ₖ ,    w̃ₖᵀ = ∑ wⱼᵀ 𝓣ⱼ,ₖ
%                ʲ⁼¹                  ʲ⁼¹
%   In the two-step approach to computing the nonlinear
%   input-normal/output-diagonal transformation, the first transformation is
%   linear, so x = Φ(z) = T₁z and the sums only contain one term each, so the
%   transformed coefficients are given by the simpler formulas
%           ṽₖ = 𝓣ₖ,ₖᵀ vₖ ,    w̃ₖ = 𝓣ₖ,ₖᵀ wₖ
%   where 𝓣ₖ,ₖ = T₁⊗T₁⊗...⊗T₁ (k times). The calTTv function helps with
%   doing this procedure efficiently using Kronecker product identities.
%
%   Authors: Nick Corbin, UCSD
%
%   See also: calTTv
%%
vec = @(X) X(:);

degree = length(w);
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
