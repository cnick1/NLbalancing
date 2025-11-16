function [Tbal] = balancingTransformation(v, w, degree, verbose)
%balancingTransformation Return a polynomial balancing transformation x = ̅Φ(z̄) = Φ(𝝋(z̄))
%
%   Usage:  Tbal = balancingTransformation(TinOd, Tscal)
%
%   Inputs:
%       v,w     - cell arrays containing the polynomial energy function
%                 coefficients; these should already be in input-normal form.
%       degree  - desired degree of the computed transformation (default =
%                 degree of energy functions - 1).
%       verbose - optional argument to print runtime information.
%
%   Outputs:     Tbal - cell array containing balancing transformation
%                       coefficients                              ( ̅Φ(z̄) )
%
%   Description: The nonlinear balancing transformation x = ̅Φ(z̄) = Φ(𝝋(z̄))
%   is the composition of the input-normal/output-diagonal transformation
%   x = Φ(z) and the scaling transformation z = 𝝋(z̄). The balancing
%   transformation puts the energy functions in the form
%           𝓔⁻( ̅Φ(z̄)) = 1/2 z̄ᵀ Σ⁻¹(z̄) z̄
%           𝓔⁺( ̅Φ(z̄)) = 1/2 z̄ᵀ  Σ(z̄)  z̄
%   where Σ(z̄) is the diagonal matrix of singular value functions ̅σᵢ(z̄ᵢ).
%   In this function, we compute an approximate polynomial expansion for
%   the composite transformation. This is done in three steps:
%       1) Compute the input-normal/output-diagonal transformation x = Φ(z)
%       2) Use the squared singular value functions (observability energy in
%       the input-normal/output-diagonal coordinates) to compute the
%       scaling transformation z = 𝝋(z̄)
%       3) Construct the balancing transformation as the composition of the
%       two transformations x = ̅Φ(z̄) = Φ(𝝋(z̄))
%
%   References: [1] J. M. A. Scherpen, “Balancing for nonlinear systems,”
%               Systems & Control Letters, vol. 21, no. 2, pp. 143–153, Aug.
%               1993, doi: 10.1016/0167-6911(93)90117-o.
%               [2] J. M. A. Scherpen and A. J. Van Der Schaft, “Normalized
%               coprime factorizations and balancing for unstable nonlinear
%               systems,” International Journal of Control, vol. 60, no. 6, pp.
%               1193–1222, Aug. 1993, doi: 10.1080/00207179408921517.
%               [3] J. M. A. Scherpen, “Balancing for nonlinear systems,”
%               University of Twente, 1994.
%               [4] K. Fujimoto and J. M. A. Scherpen, “Balanced realization and
%               model order reduction for nonlinear systems based on singular
%               value analysis,” SIAM Journal on Control and Optimization, vol.
%               48, no. 7, pp. 4591–4623, Jan. 2010, doi: 10.1137/070695332.
%               [5] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” IEEE
%               Transactions on Automatic Control, vol. 70, no. 5, pp.
%               3088–3100, May 2025, doi: 10.1109/tac.2024.3494472.
%               [6] N. A. Corbin, A. Sarkar, J. M. A. Scherpen, and B. Kramer,
%               “Scalable computation of input-normal/output-diagonal balanced
%               realization for control-affine polynomial systems,” Systems &
%               Control Letters, vol. 204, p. 106178, Oct. 2025, doi:
%               10.1016/j.sysconle.2025.106178.
%
%   Part of the NLbalancing repository.
%
%   See also: approxPastEnergy, approxFutureEnergy, inputNormalOutputDiagonalTransformation, composePolynomials
%%
arguments
    v
    w
    degree = length(v) - 1
    verbose = false
end

[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree, verbose);
[Tscal, TscalInv] = scalingTransformation(sigmaSquared, degree);
Tbal = composePolynomials(TinOd, Tscal);

end

function [Tscal, TscalInv] = scalingTransformation(sigmaSquared, d)
%scalingTransformation Return polynomial expansions for z = 𝝋(z̄) and ̄z = 𝝋⁻¹(z)
% i.e. scaling transformation and its inverse
%
%   Usage:  [Tscal, TscalInv] = scalingTransformation(sigmaSquared)
%
%   Inputs:
%        sigmaSquared - the coefficients of the squared singular value functions
%             options - name-value pair optional arguments
%
%   Outputs:    Tscal - the coefficients of the scaling transformation
%            TscalInv - the coefficients of the inverse transformation
%
%   Description: The nonlinear balancing transformation is x = ̅Φ(z̄) =
%   Φ(𝝋(z̄), which composition of the input-normal/output-diagonal
%   transformation x = Φ(z) and the scaling transformation z = 𝝋(z̄), where
%           zᵢ = 𝜑ᵢ(z̄) = z̄ᵢ / √ ̅σᵢ(z̄ᵢ)
%   However, we have a polynomial expansion for σ²ᵢ(zᵢ), i.e. the squared
%   singular value functions in terms of z, not z̄.
%
%   A previous approach for computing z = 𝝋(z̄) was therefore to compute
%   the inverse of ̄z = 𝝋⁻¹(z) via Newton iteration for a single point. In
%   this function, we instead want to get a polynomial approximation for
%   both z = 𝝋(z̄) and ̄z = 𝝋⁻¹(z). We begin with the function ̄z=𝝋⁻¹(z),
%   whose components are given by
%           z̄ᵢ = 𝜑⁻¹ᵢ(z) = zᵢ ⁴√ σ²ᵢ(zᵢ)       (= zᵢ √ σᵢ(zᵢ))
%   This can be approximated via Taylor expansion using the usual formula.
%   It is a bit messy to do by hand, so for now we will use symbolic
%   computations. It is always the same though, so we could also hard-code
%   up to some degree to speed this up; I don't think it will be an issue
%   though symbolically.
%
%   Then, we can compute the polynomial expansion for the inverse
%   transformation using series reversion [1,2].
%
%   References: [1] https://mathworld.wolfram.com/SeriesReversion.html
%               [2] M. Abramowitz and I. A. Stegun, Eds., Handbook of
%               mathematical functions, 9. Dover print.; [Nachdr. der Ausg.
%               von 1972]. in Dover books on mathematics. New York, NY:
%               Dover Publ., 2013. Page 16.
%
%   Part of the NLbalancing repository.
%%
arguments
    sigmaSquared
    d = size(sigmaSquared, 2)
end

[n,nd] = size(sigmaSquared);

c = zeros(n,9);
c(:,1:nd) = sigmaSquared;

a = [c(:,1).^(1/4), ... % this is S.^(1/2) for linear balancing
    c(:,2)./(4.*c(:,1).^(3/4)), ...
    (-3.*c(:,2).^2 + 8.*c(:,1).*c(:,3))./(32.*c(:,1).^(7/4)), ...
    (7.*c(:,2).^3 - 24.*c(:,1).*c(:,2).*c(:,3) + 32.*c(:,1).^2.*c(:,4))./(128.*c(:,1).^(11/4)), ...
    1/(2048.*c(:,1).^(15/4)).*(-77*c(:,2).^4 + 336.*c(:,1).*c(:,2).^2.*c(:,3) - 384.*c(:,1).^2.*c(:,2).*c(:,4) + 64.*c(:,1).^2.*(-3*c(:,3).^2 + 8.*c(:,1).*c(:,5))), ...
    1/( 8192.*c(:,1).^(19/4)).*(231.*c(:,2).^5-1232.*c(:,1).*c(:,2).^3.*c(:,3)+1344.*c(:,1).^2.*c(:,2).^2.*c(:,4)-192.*c(:,1).^2.*c(:,2).*(-7.*c(:,3).^2+8.*c(:,1).*c(:,5))+512.*c(:,1).^3.*(-3.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))), ...
    1/(65536.*c(:,1).^(23/4)).*(-1463.*c(:,2).^6+9240.*c(:,1).*c(:,2).^4.*c(:,3)-9856.*c(:,1).^2.*c(:,2).^3.*c(:,4)+1344.*c(:,1).^2.*c(:,2).^2.*(-11.*c(:,3).^2+8.*c(:,1).*c(:,5))-3072.*c(:,1).^3.*c(:,2).*(-7.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))+512.*c(:,1).^3.*(7.*c(:,3).^3-24.*c(:,1).*c(:,3).*c(:,5)+4.*c(:,1).*(-3.*c(:,4).^2+8.*c(:,1).*c(:,7)))), ...
    1/(262144.*c(:,1).^(27/4)).*(4807.*c(:,2).^7-35112.*c(:,1).*c(:,2).^5.*c(:,3)+36960.*c(:,1).^2.*c(:,2).^4.*c(:,4)-4928.*c(:,1).^2.*c(:,2).^3.*(-15.*c(:,3).^2+8.*c(:,1).*c(:,5))+10752.*c(:,1).^3.*c(:,2).^2.*(-11.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))-512.*c(:,1).^3.*c(:,2).*(77.*c(:,3).^3-168.*c(:,1).*c(:,3).*c(:,5)+12.*c(:,1).*(-7.*c(:,4).^2+8.*c(:,1).*c(:,7)))+2048.*c(:,1).^4.*(21.*c(:,3).^2.*c(:,4)-24.*c(:,1).*c(:,3).*c(:,6)+8.*c(:,1).*(-3.*c(:,4).*c(:,5)+4.*c(:,1).*c(:,8)))), ...
    1/(8388608.*c(:,1).^(31/4)).*(-129789.*c(:,2).^8+1076768.*c(:,1).*c(:,2).^6.*c(:,3)-1123584.*c(:,1).^2.*c(:,2).^5.*c(:,4)+147840.*c(:,1).^2.*c(:,2).^4.*(-19.*c(:,3).^2+8.*c(:,1).*c(:,5))-315392.*c(:,1).^3.*c(:,2).^3.*(-15.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))+43008.*c(:,1).^3.*c(:,2).^2.*(55.*c(:,3).^3-88.*c(:,1).*c(:,3).*c(:,5)+4.*c(:,1).*(-11.*c(:,4).^2+8.*c(:,1).*c(:,7)))-49152.*c(:,1).^4.*c(:,2).*(77.*c(:,3).^2.*c(:,4)-56.*c(:,1).*c(:,3).*c(:,6)+8.*c(:,1).*(-7.*c(:,4).*c(:,5)+4.*c(:,1).*c(:,8)))+4096.*c(:,1).^4.*(-77.*c(:,3).^4+336.*c(:,1).*c(:,3).^2.*c(:,5)-48.*c(:,1).*c(:,3).*(-7.*c(:,4).^2+8.*c(:,1).*c(:,7))+64.*c(:,1).^2.*(-3.*c(:,5).^2-6.*c(:,4).*c(:,6)+8.*c(:,1).*c(:,9))))];

A = [1./a(:,1),...
    -(a(:,2)./a(:,1).^3),...
    (2.*a(:,2).^2-a(:,1).*a(:,3))./a(:,1).^5,...
    -((5.*a(:,2).^3-5.*a(:,1).*a(:,2).*a(:,3)+a(:,1).^2.*a(:,4))./a(:,1).^7),...
    1./a(:,1).^9.*(14.*a(:,2).^4-21.*a(:,1).*a(:,2).^2.*a(:,3)+3.*a(:,1).^2.*a(:,3).^2+6.*a(:,1).^2.*a(:,2).*a(:,4)-a(:,1).^3.*a(:,5)),...
    1./a(:,1).^11.*(-42.*a(:,2).^5+84.*a(:,1).*a(:,2).^3.*a(:,3)-28.*a(:,1).^2.*a(:,2).^2.*a(:,4)+7.*a(:,1).^2.*a(:,2).*(-4.*a(:,3).^2+a(:,1).*a(:,5))+a(:,1).^3.*(7.*a(:,3).*a(:,4)-a(:,1).*a(:,6))),...
    1./a(:,1).^13.*(132.*a(:,2).^6-330.*a(:,1).*a(:,2).^4.*a(:,3)+120.*a(:,1).^2.*a(:,2).^3.*a(:,4)-36.*a(:,1).^2.*a(:,2).^2.*(-5.*a(:,3).^2+a(:,1).*a(:,5))+8.*a(:,1).^3.*a(:,2).*(-9.*a(:,3).*a(:,4)+a(:,1).*a(:,6))+a(:,1).^3.*(-12.*a(:,3).^3+8.*a(:,1).*a(:,3).*a(:,5)+a(:,1).*(4.*a(:,4).^2-a(:,1).*a(:,7)))),...
    1./a(:,1).^15.*(-429.*a(:,2).^7+1287.*a(:,1).*a(:,2).^5.*a(:,3)-495.*a(:,1).^2.*a(:,2).^4.*a(:,4)+165.*a(:,1).^2.*a(:,2).^3.*(-6.*a(:,3).^2+a(:,1).*a(:,5))-45.*a(:,1).^3.*a(:,2).^2.*(-11.*a(:,3).*a(:,4)+a(:,1).*a(:,6))+3.*a(:,1).^3.*a(:,2).*(55.*a(:,3).^3-30.*a(:,1).*a(:,3).*a(:,5)+3.*a(:,1).*(-5.*a(:,4).^2+a(:,1).*a(:,7)))+a(:,1).^4.*(-45.*a(:,3).^2.*a(:,4)+9.*a(:,1).*a(:,3).*a(:,6)+a(:,1).*(9.*a(:,4).*a(:,5)-a(:,1).*a(:,8)))),...
    1/a(:,1).^17.*(1430.*a(:,2).^8-5005.*a(:,1).*a(:,2).^6.*a(:,3)+2002.*a(:,1).^2.*a(:,2).^5.*a(:,4)-715.*a(:,1).^2.*a(:,2).^4.*(-7.*a(:,3).^2+a(:,1).*a(:,5))+220.*a(:,1).^3.*a(:,2).^3.*(-13.*a(:,3).*a(:,4)+a(:,1).*a(:,6))-55.*a(:,1).^3.*a(:,2).^2.*(26.*a(:,3).^3-12.*a(:,1).*a(:,3).*a(:,5)+a(:,1).*(-6.*a(:,4).^2+a(:,1).*a(:,7)))+10.*a(:,1).^4.*a(:,2).*(66.*a(:,3).^2.*a(:,4)-11.*a(:,1).*a(:,3).*a(:,6)+a(:,1).*(-11.*a(:,4).*a(:,5)+a(:,1).*a(:,8)))+a(:,1).^4.*(55.*a(:,3).^4-55.*a(:,1).*a(:,3).^2.*a(:,5)+5.*a(:,1).*a(:,3).*(-11.*a(:,4).^2+2.*a(:,1).*a(:,7))+a(:,1).^2.*(5.*a(:,5).^2+10.*a(:,4).*a(:,6)-a(:,1).*a(:,9))))];

Tscal = cell(1,d);
TscalInv = cell(1,d);
TscalInv{1} = invertibleMatrix(diag(a(:,1)),diag(A(:,1)));
Tscal{1} = invertibleMatrix(diag(A(:,1)),diag(a(:,1)));
for i=2:d
    TscalInv{i} = reshape(sparse(linspace(1,n^(i+1),n),1,a(:,i)),n,[]);
    Tscal{i}    = reshape(sparse(linspace(1,n^(i+1),n),1,A(:,i)),n,[]);
end

end



% function coeff = phiInvCoeffs(c, i)
% % alternative way to possibly only compute up to the degree we need
% switch i
%     case 1
%         coeff = c(:,1).^(1/4);
%     case 2
%         coeff = c(:,2)./(4.*c(:,1).^(3/4));
%     case 3
%         coeff = (-3.*c(:,2).^2 + 8.*c(:,1).*c(:,3))./(32.*c(:,1).^(7/4));
%     case 4
%         coeff = (7.*c(:,2).^3 - 24.*c(:,1).*c(:,2).*c(:,3) + 32.*c(:,1).^2.*c(:,4))./(128.*c(:,1).^(11/4));
%     case 5
%         coeff = 1/(2048.*c(:,1).^(15/4)).*(-77*c(:,2).^4 + 336.*c(:,1).*c(:,2).^2.*c(:,3) - 384.*c(:,1).^2.*c(:,2).*c(:,4) + 64.*c(:,1).^2.*(-3*c(:,3).^2 + 8.*c(:,1).*c(:,5)));
%     case 6
%         coeff = 1/( 8192.*c(:,1).^(19/4)).*(231.*c(:,2).^5-1232.*c(:,1).*c(:,2).^3.*c(:,3)+1344.*c(:,1).^2.*c(:,2).^2.*c(:,4)-192.*c(:,1).^2.*c(:,2).*(-7.*c(:,3).^2+8.*c(:,1).*c(:,5))+512.*c(:,1).^3.*(-3.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6)));
%     case 7
%         coeff = 1/(65536.*c(:,1).^(23/4)).*(-1463.*c(:,2).^6+9240.*c(:,1).*c(:,2).^4.*c(:,3)-9856.*c(:,1).^2.*c(:,2).^3.*c(:,4)+1344.*c(:,1).^2.*c(:,2).^2.*(-11.*c(:,3).^2+8.*c(:,1).*c(:,5))-3072.*c(:,1).^3.*c(:,2).*(-7.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))+512.*c(:,1).^3.*(7.*c(:,3).^3-24.*c(:,1).*c(:,3).*c(:,5)+4.*c(:,1).*(-3.*c(:,4).^2+8.*c(:,1).*c(:,7))));
%     case 8
%         coeff = 1/(262144.*c(:,1).^(27/4)).*(4807.*c(:,2).^7-35112.*c(:,1).*c(:,2).^5.*c(:,3)+36960.*c(:,1).^2.*c(:,2).^4.*c(:,4)-4928.*c(:,1).^2.*c(:,2).^3.*(-15.*c(:,3).^2+8.*c(:,1).*c(:,5))+10752.*c(:,1).^3.*c(:,2).^2.*(-11.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))-512.*c(:,1).^3.*c(:,2).*(77.*c(:,3).^3-168.*c(:,1).*c(:,3).*c(:,5)+12.*c(:,1).*(-7.*c(:,4).^2+8.*c(:,1).*c(:,7)))+2048.*c(:,1).^4.*(21.*c(:,3).^2.*c(:,4)-24.*c(:,1).*c(:,3).*c(:,6)+8.*c(:,1).*(-3.*c(:,4).*c(:,5)+4.*c(:,1).*c(:,8))));
%     case 9
%         coeff = 1/(8388608.*c(:,1).^(31/4)).*(-129789.*c(:,2).^8+1076768.*c(:,1).*c(:,2).^6.*c(:,3)-1123584.*c(:,1).^2.*c(:,2).^5.*c(:,4)+147840.*c(:,1).^2.*c(:,2).^4.*(-19.*c(:,3).^2+8.*c(:,1).*c(:,5))-315392.*c(:,1).^3.*c(:,2).^3.*(-15.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))+43008.*c(:,1).^3.*c(:,2).^2.*(55.*c(:,3).^3-88.*c(:,1).*c(:,3).*c(:,5)+4.*c(:,1).*(-11.*c(:,4).^2+8.*c(:,1).*c(:,7)))-49152.*c(:,1).^4.*c(:,2).*(77.*c(:,3).^2.*c(:,4)-56.*c(:,1).*c(:,3).*c(:,6)+8.*c(:,1).*(-7.*c(:,4).*c(:,5)+4.*c(:,1).*c(:,8)))+4096.*c(:,1).^4.*(-77.*c(:,3).^4+336.*c(:,1).*c(:,3).^2.*c(:,5)-48.*c(:,1).*c(:,3).*(-7.*c(:,4).^2+8.*c(:,1).*c(:,7))+64.*c(:,1).^2.*(-3.*c(:,5).^2-6.*c(:,4).*c(:,6)+8.*c(:,1).*c(:,9))));
% end
% end

% function coefficients = varphiInv(sigmaSquared, d)
% % alternative way to do things symbolically
% %%
% arguments
%     sigmaSquared
%     d = size(sigmaSquared, 2)
% end
%
% [n,nd] = size(sigmaSquared);
% coefficients = zeros(n, d+1);
%
% syms z c0
% c = [c0, sym('c', [1, nd-1])];
%
% symvarphiInv = z * sum(c .* z.^(0:nd-1))^(1/4);
% expansion = taylor(symvarphiInv, 'Order', d);
%
% % Option 1: not vectorized
% % for i=1:n
% %     temp = double(coeffs(subs(expansion,c,sigmaSquared(i,:)),z,'All'));
% %     coefficients(i,:) = flip([zeros(1,d+1-length(temp)), temp]);
% % end
%
% % Option 2: vectorized but possibly buggy
% coefficients = cell2mat(arrayfun(@(i) coeffs2(subs(expansion,c,sigmaSquared(i,:)), z, d+1), 1:n, 'UniformOutput', false).');
%
% coefficients = coefficients(:,2:end);
%
% end
%
% function c = coeffs2(p,var,d)
% c = zeros(1,d);
% temp = coeffs(p,var,'All');
% c(1:length(temp)) = double(flip(temp));
% end