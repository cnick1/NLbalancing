function [F] = composePolynomials(P,T,varargin,nvp)
%composePolynomials Combines polynomial into one: f(z) = p(Φ(...Tₙ(z)))
%
%   Usage: [F] = composePolynomials(P,T)
%
%   Inputs:
%       P - cell array containing the polynomial coefficients for the first
%           function p(x)
%       T - cell array containing the polynomial coefficients for the second
%           function Φ(z)
%       varargin - the function permits passing 2 or more transformations; if
%                  more than two are passed, the function will compute the
%                  result using recursion.
%     Optional name/value pair inputs:
%      degree - desired degree of output. By default, we will compute the full
%               result, which is the degree of p(x) multiplied by the degree of Φ(z)
%
%   Output:
%       F - cell array containing the polynomial coefficients for the
%           composed function.
%
%   Note: Generally, combining transformations increases the polynomial degree
%   of the resulting transformation, so this is not always the most efficient
%   way to do things. For example, if you combine two degree 3 transformations,
%   you get a degree 9 transformation. It would likely be more efficient to just
%   evaluate them one after the other.
%
%   Description: (Lemma 1, [1]) Consider the function p(Φ(z)), where
%           p(x) = P₁x + P₂(x⊗x) + ... + Pd(x...⊗x)               (1)
%           Φ(z) = T₁z + T₂(z⊗z) + ... + Tₖ(z...⊗z)                (2)
%   In the new z coordinates, f(z) = p(Φ(z)) can be written to degree d as
%           f(z) = F₁z + F₂(z⊗z) + ... + Fd(z...⊗z)
%   where             ____________________
%                    |        ᵢ           |
%                ->  |   Fᵢ = ∑ Pⱼ 𝓣ⱼ,ᵢ   |                        (3)
%                    |_______ʲ⁼¹__________|
%
%   and 𝓣ⱼ,ᵢ is calligraphic T, used to compactly write the sum of all unique
%   tensor products with j factors and nⁱ columns [2]. This is easily
%   verifiable: inserting the expansions (1) and (2) gives
%           p(Φ(z)) = P₁(Φ(z)) + P₂(Φ(z)⊗Φ(z)) + ...
%                   = P₁(T₁z + T₂(z⊗z) + ...)
%                     + P₂( (T₁z + T₂(z⊗z) + ...) ⊗ (T₁z + T₂(z⊗z) + ...) )
%                     + ...
%   This is exactly what is implemented in this function. in the derivation
%   above.  The quantities Pⱼ 𝓣ⱼ,ᵢ can be computed using the function calTTv()
%   from the KroneckerTools repository; it is possible that improvements can be
%   made to avoid having to transpose things.
%
%   References: [1] N. A. Corbin, A. Sarkar, J. M. A. Scherpen, and B. Kramer,
%                “Scalable computation of input-normal/output-diagonal balanced
%                realization for control-affine polynomial systems,” Systems &
%                Control Letters, vol. 204, p. 106178, Oct. 2025, doi:
%                10.1016/j.sysconle.2025.106178.
%               [2] B. Kramer, S. Gugercin, and J. Borggaard, “Nonlinear
%                balanced truncation: Part 2—model reduction on manifolds,”
%                arXiv, Feb. 2023. doi: 10.48550/ARXIV.2302.02036
%
%   Authors: Nick Corbin, UCSD
%
arguments
    P cell
    T cell
end
arguments (Repeating)
    varargin
end
arguments
    % nvp.degree = length(P)
    nvp.degree = length(P) * length(T)
end

[n1,n2] = size(P{1});
[n2,n3] = size(T{1});

F = cell(1,nvp.degree);

F{1} = calTTv(T, 1, 1, P{1}.').';
for i = 2:nvp.degree
    F{i} = sparse(n1,n3^i);
    % This next loop evaluates the summation in equation (3)
    for j=flip(1:i)
        if j > length(P); continue; end   % skip if requesting P that doesn't exist
        if (i-j)-length(T)>=0; break; end % break when we run out of Ts
        F{i} = F{i} + calTTv(T, j, i, P{j}.').';
    end
end

if nargin > 2 % Compute the remaining ones by recursion
    [F] = composePolynomials(F,varargin{3},varargin{4:end},degree=nvp.degree);
end

end
