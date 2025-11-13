function [Tout] = composeTransformations(T1,T2,varargin)
%composeTransformations Combines polynomial transformations into one: x = Tᵒᵘᵗ(z) = T₁(T₂(...Tₙ(z)))
% DEPRECATED; replace with composePolynomials instead
%
%   Usage: [Tout] = composeTransformations(T1,T2)
%
%   Inputs:
%       T1,T2,... - cell arrays containing the polynomial transformation
%                   coefficients. At least two transformations are required,
%                   but more can be included and the function will apply
%                   the transformations one at a time recursively.
%
%   Output:
%       Tout - cell array containing the polynomial coefficients for the
%              composed transformation.
%
%   Note: Generally, combining transformations increases the polynomial degree
%   of the resulting transformation, so this is not always the most efficient
%   way to do things. For example, if you combine two degree 3 transformations,
%   you get a degree 9 transformation. It would likely be more efficient to just
%   evaluate them one after the other. Additionally, in this function, we
%   truncate the result to the degree of the first transformation, so there will
%   be a truncation error.
%
%   Description: (Lemma 1, [1]) Consider the function p(Φ(z)), where
%           p(x) = P₁x + P₂(x⊗x) + ... + Pd(x...⊗x)                      (1)
%           Φ(z) = T₁z + T₂(z⊗z) + ... + Tₖ(z...⊗z)                      (2)
%   In the new z coordinates, f(z) = p(Φ(z)) can be written to degree d as
%           f(z) = F₁z + F₂(z⊗z) + ... + Fd(z...⊗z)
%   where
%            ᵢ
%       Fᵢ = ∑ Pⱼ 𝓣ⱼ,ᵢ                                                    (3)
%           ʲ⁼¹
%   and 𝓣ⱼ,ᵢ is calligraphic T, used to compactly write the sum of all unique
%   tensor products with j factors and nⁱ columns [2]. This is easily
%   verifiable: inserting the expansions (1) and (2) gives
%           p(Φ(z)) = P₁(Φ(z)) + P₂(Φ(z)⊗Φ(z)) + ...
%                   = P₁(T₁z + T₂(z⊗z) + ...)
%                     + P₂( (T₁z + T₂(z⊗z) + ...) ⊗ (T₁z + T₂(z⊗z) + ...) )
%                     + ...
%   This is exactly what is implemented in this function, where p(x) = T₁(x) and
%   Φ(z) = T₂(z) in the derivation above. See also transformDynamics() for
%   possible improvements. The quantities Pⱼ 𝓣ⱼ,ᵢ can be computed using the
%   function calTTv() from the KroneckerTools repository; it is possible that
%   improvements can be made to avoid having to transpose things. This function
%   could also be adapted to take an optional polynomial degree argument.
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
error('This function is deprecated, use <a href="matlab:open composePolynomials">composePolynomials</a> instead.')
ld = length(T1);
n = length(T1{1});

Tout = cell(size(T1));


Tout{1}= calTTv(T2, 1, 1, T1{1}.').';
for i = 2:ld
    Tout{i} = zeros(n,n^i);
    for j=1:i  % This is the summation in equation (3)
        Tout{i}= Tout{i} + calTTv(T2, j, i, T1{j}.').';
    end
end

if nargin > 2 % Compute the remaining ones by recursion
    [Tout] = composeTransformations(Tout,varargin{3},varargin{4:end});
end

end
