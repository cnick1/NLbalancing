function P = transformationInverse(T)
%transformationInverse Return the coefficients of the inverse of the transformation.
%
%   Usage:  P = transformationInverse(T)
%
%   Inputs:      T - cell array containing transformation coefficients
%
%   Outputs:     P - cell array of the inverse transformation coefficients
%
%   Description: Consider a transformation x = Φ(z), for which we have a
%   polynomial approximation of the form
%
%            Φ(z)  =  T₁z + T₂(z⊗z) + ... + Tₖ(z...⊗z)                   (1)
%
%   We wish to find a polynomial approximation to the inverse transformation
%   z = Φ⁻¹(x), which we express as
%
%          Φ⁻¹(x)  =  P₁x + P₂(x⊗x) + ... + Pₖ(x...⊗x)                   (2)
%
%   By the definition of the inverse, Φ⁻¹(Φ(z)) = z, so by following a similar
%   procedure to that described in composePolynomials, we can derive the
%   expressions for the inverse transformation coefficients P₁, P₂,..., Pₖ. The
%   result of equating sets of terms of the same degree is
%         O(z):   I = P₁𝓣₁,₁ = P₁T₁       ->   P₁ = T₁⁻¹
%                   :
%        O(z³):   0 = P₁𝓣₁,₃ + P₂𝓣₂,₃ + P₃𝓣₃,₃
%                   :
%                   : ᵢ           ᵢ₋₁
%        O(zⁱ):   0 = ∑ Pⱼ 𝓣ⱼ,ᵢ = ∑ Pⱼ 𝓣ⱼ,ᵢ + Pᵢ 𝓣ᵢ,ᵢ
%                    ʲ⁼¹          ʲ⁼¹
%              _______________________________
%             |         ᵢ₋₁                   |
%         ->  |  Pᵢ = -( ∑ Pⱼ 𝓣ⱼ,ᵢ ) 𝓣ᵢ,ᵢ⁻¹  |                          (3)
%             |_________ʲ⁼¹___________________|
%
%   In general, Pᵢ ≠ Tᵢ⁻¹, so 𝓣ₘ,ₖ⁻¹ ≠ 𝓟ₘ,ₖ; however, since P₁ = T₁⁻¹, it is
%   true that 𝓣ᵢ,ᵢ⁻¹ = 𝓟ᵢ,ᵢ. This is taken advantage of in the last step to
%   replace the inversion on the right with regular matrix multiplication.
%
%   A critical detail here: the only inversions required are of T₁, which is
%   assumed to be invertible. Often, such as the case of balanced truncation,
%   the inverse T₁⁻¹ may have a closed-form analytical representation, so we
%   don't need to actually compute it or solve any linear systems; everything
%   can be done with matrix multiplication.
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
P = cell(size(T)); % default to producing an expansion the same size as T

n = size(T{1}, 1);

P{1} = invertibleMatrix(inv(T{1}), T{1}); % See invertibleMatrix; if T{1} is an invertibleMatrix already, inv() does not actually have to compute it
for i=2:length(T)
    P{i} = zeros(size(T{i}));
    for j = 1:(i-1) % Compute the sum ∑ Pⱼ 𝓣ⱼ,ᵢ
        % for idx = 1:n
        %     P{i}(idx,:) = P{i}(idx,:) - calTTv(T,j,i,P{j}(idx,:).').';
        % end
        P{i} = P{i} - calTTv(T,j,i,P{j}.').';
    end
    % Now multiply by the inverse on right 𝓣ᵢ,ᵢ⁻¹
    % for idx = 1:n
    %     P{i}(idx,:) = calTTv(P,i,i,P{i}(idx,:).').';
    % end
    P{i} = calTTv(P,i,i,P{i}.').';
    
    % Row symmetrization may be necessary
end

end