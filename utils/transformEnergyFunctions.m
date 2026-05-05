function [vtilde, wtilde] = transformEnergyFunctions(v, w, T, inputNormal)
%transformEnergyFunctions Transforms the energy coefficients v and w by T.
%
%   Usage: [vtilde, wtilde] = transformEnergyFunctions(v, w, T)
%
%   Inputs:
%       v,w         - cell arrays containing the polynomial energy function
%                     coefficients
%       T           - cell array containing the polynomial transformation
%                     coefficients
%       inputNormal - boolean indicating the transformation is input
%                     normal, used to speed up the computation by skipping
%                     that transformation (default = false)
%
%   Output:
%       v,w   - cell arrays containing the transformed polynomial energy function
%               coefficients
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
%
%   Authors: Nick Corbin, UCSD
%
if nargin < 4
    inputNormal = false;
end

vec = @(X) X(:);

degree = length(w);
% n = sqrt(length(v{2}));
[n, r] = size(T{1});
V2 = reshape(v{2}, n, n);
W2 = reshape(w{2}, n, n);

vtilde = cell(1, degree);
wtilde = cell(1, degree);

vtilde{2} = vec(T{1}.' * V2 * T{1});
wtilde{2} = vec(T{1}.' * W2 * T{1});

if inputNormal
    % TODO: consider adding a check to verify that the transformation is input normal, throw error if not
    for k = 3:degree
        vtilde{k} = sparse(n^k,1);
        try
            wtilde{k} = vec(T{k - 1}.' * W2 * T{1}) + vec(T{1}.' * W2 * T{k - 1});
        catch
        end

        for i = 2:k - 2
            j = k - i;
            try
                wtilde{k} = wtilde{k} + vec(T{j}.' * W2 * T{i});
            catch
            end
        end

        for i = 3:k
            try
                wtilde{k} = wtilde{k} + calTTv(T, i, k, w{i});
            catch
            end
        end

        wtilde{k} = kronMonomialSymmetrize(wtilde{k},n,k);
    end
else
    for k = 3:degree
        vtilde{k} = zeros(size(v{k}));
        wtilde{k} = zeros(size(w{k}));
        try % cursed, need to fix
            vtilde{k} = vec(T{k - 1}.' * V2 * T{1}) + vec(T{1}.' * V2 * T{k - 1});
            wtilde{k} = vec(T{k - 1}.' * W2 * T{1}) + vec(T{1}.' * W2 * T{k - 1});
        catch
        end

        for i = 2:k - 2
            j = k - i;
            try
                vtilde{k} = vtilde{k} + vec(T{j}.' * V2 * T{i});
                wtilde{k} = wtilde{k} + vec(T{j}.' * W2 * T{i});
            catch
            end
        end

        for i = 3:k
            try
                vtilde{k} = vtilde{k} + calTTv(T, i, k, v{i});
                wtilde{k} = wtilde{k} + calTTv(T, i, k, w{i});
            catch
            end
        end

        vtilde{k} = kronMonomialSymmetrize(vtilde{k},r,k);
        wtilde{k} = kronMonomialSymmetrize(wtilde{k},r,k);
    end
end
end
