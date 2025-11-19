function [ft, gt, ht] = transformDynamics(f, g, h, T, nvp)
%transformDynamics Transform the control-affine dynamic model given by f, g, h by transformation T.
%   This function returns the expansions for the transformed dynamics in the form
%       ż = f̃(z) + g̃(z) u, y = h̃(z)
%
%   Usage: [ft, gt, ht] = transformDynamics(f, g, h, T)
%
%   Inputs:
%       f,g,h - cell arrays containing the polynomial coefficients for
%               the drift, input, and output in the original coordinates.
%       T     - cell array containing the polynomial transformation
%               coefficients.
%
%   Output:
%       ft,gt,ht - cell arrays containing the polynomial coefficients
%                  for the transformed drift, input, and output.
%
%   Description: Given a transformation x = Φ(z), we seek to represent the
%   dynamics for the control-affine system in the new coordinates as
%            ẋ = f(x) + g(x) u,      ->       ż = f̃(z) + g̃(z) u,
%            y = h(x),               ->       y = h̃(z).
%   In general, it may not be possible to do this exactly. Applying the
%   transformation yields
%           ∂Φ(z)/∂z ż = f(Φ(z)) + g(Φ(z)) u                         (1)
%                    y = h(Φ(z))
%   and we do not in general have an exact way to write [∂Φ(z)/∂z]⁻¹. In this
%   function, we will approximate the functions f̃(z), g̃(z), h̃(z) by computing
%   their Taylor expansions. This can be done exactly for linear transformations
%   but for polynomial approximations involves a truncation of the expansions.
%   Instead of computing [∂Φ(z)/∂z]⁻¹ explicitly, we will leave it on the
%   left-hand-side and match terms of the same degree; doing it this way, we
%   only ever have to invert T₁, which we may have a closed-form analytical
%   expression for.
%
%   Inserting ż = f̃(z) + g̃(z) u into (1), we obtain
%           ∂Φ(z)/∂z (f̃(z) + g̃(z) u) = f(Φ(z)) + g(Φ(z)) u           (2)
%                                  y = h(Φ(z))
%   From here, we can expand everything as polynomials and match terms of
%   the same degree:
%           ∂Φ(z)/∂z f̃(z) = f(Φ(z)),
%           ∂Φ(z)/∂z g̃(z) = g(Φ(z)),
%                    h̃(z) = h(Φ(z)).
%   The polynomial expansions for these quantities are 
%           ∂Φ(z)/∂z = T₁ + 2T₂(I⊗x) + ... + d Td(I...⊗x)
%           f̃(z)     = F̃₁x + F̃₂(z⊗z) + ... + F̃d(z...⊗z)
%           f(Φ(z))  = F₁z + F₂(z⊗z) + ... + Fd(z...⊗z)
%           g̃(z)     = G̃₁x + G̃₂(z⊗z) + ... + G̃d(z...⊗z)
%           g(Φ(z))  = G₁z + G₂(z⊗z) + ... + Gd(z...⊗z)
%           h̃(z) = h(Φ(z)) = H₁z + H₂(z⊗z) + ... + Hd(z...⊗z)
%
%   After matching the terms of the same degree, the following expressions
%   are found for the coefficients F̃ₖ and G̃ₖ:
%                __________________________________
%               |                ₖ                 |
%               |   F̃ₖ = T⁻¹(Fₖ - ∑ Tᵢℒᵢ(F̃ₖ₋ᵢ₋₁)   |                  (3)
%         ->    |_______________ⁱ⁼²_______________|
%               |               ₖ₊₁                |
%               |   G̃ₖ = T⁻¹(Gₖ - ∑ Tᵢℒᵢ(G̃ₖ₋ᵢ₋₁)   |                  (4)
%               |_______________ⁱ⁼²_______________|
%
%   The procedure is therefore:
%       Step 1) Compute f(Φ(z)), g(Φ(z)), h(Φ(z)) as the compositions
%               of polynomials; after this, we already have h̃(z)
%       Step 2) Compute the coefficients of f̃(z), g̃(z) by matching
%               terms on the left and right
%
%   Authors: Nick Corbin, UCSD
%
%  See also: invertibleMatrix
%%
arguments
    f
    g
    h
    T
    nvp.degree = length(f)
    nvp.r = size(T{1},2)
end
vec = @(X) X(:);

ld = length(T);
[n,m] = size(g{1});

lf = length(f); lg = length(g);

%% Step 1) Compute transformed f(Φ(z)), g(Φ(z)), h̃(z) = h(Φ(z))
% Transform drift f(x)
ft = composePolynomials(f,T,degree=nvp.degree);

% Transform input map g(x)
gt = cell(1,max(max(1,lg-1)*ld,lf*ld));
gt{1} = g{1}; % B is not state dependent

if lg > 1 % Else if g(x) = B, we're done
    % Convert g(x) to ∑ gᵢ(x)
    g_i = cell(m,lg*ld);
    for i=1:m
        for j=1:(lg-1)
            g_i{i,j+1} = g{j+1}(:,i:m:end);
        end
    end

    % Transform gᵢ(x)
    for i=1:m
        gtemp = composePolynomials(g_i(i,2:lg),T,degree=nvp.degree);
        g_i(i,2:length(gtemp)+1) = gtemp;
    end

    % Now convert transformed ∑ gᵢ(x) vectors back to g(x) matrix
    for j=1:min((lg-1),nvp.degree)
        gt{j+1} = zeros(n,m*nvp.r^j); % The final answer is r x mr^j, but that happens after multiplying the Jacobian inverse 
        for i=1:m
            gt{j+1}(:,i:m:end) = g_i{i,j+1};
        end
    end
end

% Transform output map h(x)
ht = composePolynomials(h,T,degree=nvp.degree);

%% Step 2) Compute f̃(z), g̃(z)
% Rather than multiply by the inverse of the Jacobian, we do a coefficient
% matching that leads to the following sets of terms that need to be combined.
%  Note: If the analytical inverse of T{1} is known, use the invertibleMatrix
%  class to pass it, which overloads the matrix inversion

% Decide on the order of the dynamics
if isempty(nvp.degree)
    nvp.degree = lf; % The transformed ft, gt above are already this high
    % nvp.degree = max(lf*ld, (lg-1)*ld); % The transformed ft, gt above are already this high
    % nvp.degree = max(lf*ld*ld, (lg-1)*ld*ld); % It is known that the degree can be at least this high
    if nvp.degree > 1
        warning('Polynomial truncation degree set automatically; there may be neglected higher-order terms. ')
    end
end

% Drift f(x)
for k = 1:nvp.degree
    if k > length(ft) || isempty(ft{k})
        ft{k} = sparse(n,nvp.r^k);
    end
    for i = flip(2:k)                                       % Theoretical sum limits for F_p's; index backwards in i so that p indexes forwards
        p = k - i + 1;
        if i > ld; continue; end                            % Only run if Ti exists
        if p > length(ft); break; end                       % Only run while we have F_p's left
        %%% Naive
        for ii=1:n
            ft{k}(ii,:) = ft{k}(ii,:) - LyapProduct(ft{p}.', T{i}(ii,:).', i).';
        end
    end
    ft{k} = T{1}\ft{k};
end

% Input map g(x)
for k = 0:nvp.degree-1                                      % lg-1 to deal with zero indexing of Gs
    if k+1 > length(gt) || isempty(gt{k+1})
        gt{k+1} = sparse(n,m*nvp.r^k);
    end
    for i = flip(2:k+1)                                     % Theoretical sum limits for G_p's; index backwards in i so that p indexes forwards
        for jj=1:m
            p = k - i + 1 +1;                               % Additional +1 to deal with zero indexing of Gs
            if i > ld; continue; end                        % Only run if Ti exists
            if p > length(gt); break; end                   % Only run while we have G_p's left
            %%% Naive
            for ii=1:n
                gt{k+1}(ii,jj:m:end) = gt{k+1}(ii,jj:m:end) - LyapProduct(gt{p}(:,jj:m:end).', T{i}(ii,:).', i).';
            end
        end
    end
    gt{k+1} = T{1}\gt{k+1};
end

% output equation does not get this treatment

%% Truncate to the desired degree
% Could probably avoid computing anything higher than the desired degree if
% it is provided
ft = ft(1:min(nvp.degree,length(ft)));
gt = gt(1:min(nvp.degree,length(gt)));
ht = ht(1:min(nvp.degree,length(ht)));
end
