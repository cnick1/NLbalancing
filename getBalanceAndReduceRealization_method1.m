function [fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceAndReduceRealization(f,g,h,nvp)
%getBalanceAndReduceRealization Returns a balanced ROM for f,g,h in
%Kronecker polynomial form using simultaneous balance AND reduce. 
% This approach fails: it appears that the components of the linear
% transformation which are truncated in the input-normal/output-diagonal
% transformation are necessary for accurately computing the higher-order
% terms, so although it is mathematically possible to truncate T1 and then
% compute the higher-order terms, it is not the same as computing the full
% transformation and then truncating.
%
%   Usage: [fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceAndReduceRealization(f,g,h)
%
%   Inputs:
%       f,g,h   - cell arrays containing the polynomial coefficients
%                 for the drift, input, and output.
%                   • f must contain at least a linear drift  (A matrix)
%                   • g must contain at least a linear input  (B matrix)
%                   • h must contain at least a linear output (C matrix)
%
%     Optional name/value pair inputs:
%           eta - η=1-1/γ², where γ is the H∞ gain parameter.
%                  • For open-loop balancing, use eta=0.
%                  • For closed-loop (HJB) balancing, use eta=1.
%                 Any other value between -1 and ∞ corresponds to H∞
%                 balancing. The default is 0 for open-loop balancing.
%        degree - desired degree of the balanced realization. Ex. for linear
%                 dynamics, choosing degree=1 will produce quadratic
%                 approximations for the energy functions, linear approximations
%                 for the transformations, and the output will be a linear
%                 balanced realization. The default will be the degree of the
%                 drift term f(x) in the dynamics, i.e. we will balance
%                 everything available, no more and no less.
%        r      - reduced-order dimension if computing a balance-and-reduce
%                 transformation.
%        transformationDegree - desired degree of balancing transformation.
%                 Generally, this should be determined by the degree of
%                 balanced realization desired, i.e. the degree of the
%                 transformation through the option "degree". However, this
%                 additional option permits the user to override and for
%                 example use a linear transformation while producing a
%                 nonlinear transformed model. The most likely use case
%                 therefore corresponds to using transformationDegree=1 to
%                 compute quadratic energy functions and a linear
%                 transformation.
%       verbose - optional argument to print runtime information
%
%   Output:
%       fbal,gbal,hbal - cell arrays containing the polynomial coefficients for
%                        the drift, input, and output.
%       Tbal           - cell array containing the polynomial coefficients
%                        for the balancing transformation.
%       sigmaSquared   - squared singular value functions
%
%   Description: For a control affine dynamical system, we wish to compute a
%   balanced realization also in control-affine form via the balancing
%   transformation  x = ̅Φ(z̄):
%            ẋ = f(x) + g(x) u,      ->       ż̄ = f̃(z̄) + g̃(z̄) u
%            y = h(x).               ->       y = h̃(z̄)
%   The balanced realization is computed in three steps using helper functions:
%       1) Compute the balancing energy functions
%       2) Compute the balancing transformation
%       3) Compute the realization of the dynamics in the balanced coordinates
%   Additional details on each step can be found in the respective sections,
%   with further details in the helper functions.
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
%   See also: approxPastEnergy, approxFutureEnergy, balancingTransformation, transformDynamics
arguments
    f
    g
    h
    nvp.degree = length(f)
    nvp.r = length(f{1});
    nvp.transformationDegree = []
    nvp.eta = 0
    nvp.verbose = false
end
if isempty(nvp.transformationDegree)
    nvp.transformationDegree = nvp.degree;
end
n = length(f{1});

if nvp.verbose
    fprintf("    Drift Dynamics:\n      ")
    dispKronPoly(f)
end

%% Step 1) Compute the energy functions
% The H∞ past and future energy functions (generalizations of controllability
% and observability) are defined as
%           𝓔⁻(x) = minᵤ J(x,u) = ½∫ η||y||² + ||u||²   dt
%           𝓔⁺(x) = minᵤ J(x,u) = ½∫  ||y||² + ||u||²/η dt
% where -1 < η < ∞.  Open-loop balancing is achieved with η = 0, and closed-loop
% balancing is achieved with η = 1. The functions approxPastEnergy() and
% approxFutureEnergy() compute polynomial approximations to the energy functions
% in Kronecker product form.
[v] = approxPastEnergy(f, g, h, eta=nvp.eta, degree=nvp.transformationDegree+1, verbose=nvp.verbose);
[w] = approxFutureEnergy(f, g, h, eta=nvp.eta, degree=nvp.transformationDegree+1, verbose=nvp.verbose);

if nvp.verbose
    fprintf("    Energy functions:\n      ")
    dispKronPoly(v,n=n), fprintf("\b")
    dispKronPoly(w,n=n)
end

%% Step 2) Compute the balancing transformation
% The nonlinear balancing transformation x = ̅Φ(z̄) = Φ(𝝋(z̄)) is the
% composition of the input-normal/output-diagonal transformation x = Φ(z) and
% the scaling transformation z = 𝝋(z̄). The balancing transformation puts the
% energy functions in the diagonal balanced form
%           𝓔⁻( ̅Φ(z̄)) = 1/2 z̄ᵀ Σ⁻¹(z̄) z̄
%           𝓔⁺( ̅Φ(z̄)) = 1/2 z̄ᵀ  Σ(z̄)  z̄
% where Σ(z̄) is the diagonal matrix of singular value functions ̅σᵢ(z̄ᵢ). The
% function balancingTransformation() computes an approximate polynomial
% expansion for the balancing transformation.
[Tbal, sigmaSquared] = balancingTransformationReduced(v, w, degree=nvp.degree, verbose=nvp.verbose, r=nvp.r);

if nvp.verbose
    fprintf("    Balancing transformation:\n      ")
    dispKronPoly(Tbal)
end

%% Step 3) Compute the transformed dynamics
% Given the balancing transformation x = ̅Φ(z̄, now we seek to represent the
% dynamics for the control-affine system
%        ẋ = f(x) + g(x) u, y = h(x)
% in the new coordinates as
%        ż̄ = f̃(z̄) + g̃(z̄) u, y = h̃(z̄)
% The function transformDynamics() computes an approximate polynomial expansion
% for the transformed dynamics without inverting the nonlinear Jacobian.
[fbal,gbal,hbal] = transformDynamics(f,g,h,Tbal,degree=nvp.degree);

if nvp.verbose
    fprintf("    Transformed Drift Dynamics:\n      ")
    dispKronPoly(fbal)
end

end



function [Tbal, sigmaSquared] = balancingTransformationReduced(v, w, nvp)
%balancingTransformationReduced Return a polynomial balancing transformation x = ̅Φ(z̄) = Φ(𝝋(z̄))
%
%   Usage:  Tbal = balancingTransformationReduced(v, w)
%
%   Inputs:
%       v,w     - cell arrays containing the polynomial energy function coefficients.
%
%     Optional name/value pair inputs:
%        degree - desired degree of the balanced realization. Ex. for linear
%                 dynamics, choosing degree=1 will produce quadratic
%                 approximations for the energy functions, linear approximations
%                 for the transformations, and the output will be a linear
%                 balanced realization. The default will be on less than the
%                 degree of the v, i.e. we will balance everything available, no
%                 more and no less.
%       r       - reduced-order dimension if computing a balance-and-reduce
%                 transformation.
%       verbose - optional argument to print runtime information
%
%   Outputs:     Tbal - cell array containing balancing transformation
%                       coefficients                              ( ̅Φ(z̄) )
%        sigmaSquared - squared singular value functions
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
    v cell
    w cell
    nvp.degree = length(v) - 1
    nvp.r = sqrt(numel(v{2}))
    nvp.verbose = false
end

[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree=nvp.degree, verbose=nvp.verbose, r=nvp.r);
[Tscal, TscalInv] = scalingTransformation(sigmaSquared, degree=nvp.degree);
Tbal = composePolynomials(TinOd, Tscal, degree=nvp.degree);

end

function [Tscal, TscalInv] = scalingTransformation(sigmaSquared, nvp)
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
    nvp.degree = size(sigmaSquared, 2)
end

[n,nd] = size(sigmaSquared);

c = zeros(n,9);
c(:,1:nd) = sigmaSquared;
c1inv = diag(pinv(diag(c(:,1))));

a = [c(:,1).^(1/4), ... % this is S.^(1/2) for linear balancing
    c(:,2).*c1inv.^(3/4)./(4), ...
    (-3.*c(:,2).^2 + 8.*c(:,1).*c(:,3)).*c1inv.^(7/4)./(32), ...
    (7.*c(:,2).^3 - 24.*c(:,1).*c(:,2).*c(:,3) + 32.*c(:,1).^2.*c(:,4)).*c1inv.^(11/4)./(128), ...
    c1inv.^(15/4)/(2048).*(-77*c(:,2).^4 + 336.*c(:,1).*c(:,2).^2.*c(:,3) - 384.*c(:,1).^2.*c(:,2).*c(:,4) + 64.*c(:,1).^2.*(-3*c(:,3).^2 + 8.*c(:,1).*c(:,5))), ...
    c1inv.^(19/4)/( 8192).*(231.*c(:,2).^5-1232.*c(:,1).*c(:,2).^3.*c(:,3)+1344.*c(:,1).^2.*c(:,2).^2.*c(:,4)-192.*c(:,1).^2.*c(:,2).*(-7.*c(:,3).^2+8.*c(:,1).*c(:,5))+512.*c(:,1).^3.*(-3.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))), ...
    c1inv.^(23/4)/(65536).*(-1463.*c(:,2).^6+9240.*c(:,1).*c(:,2).^4.*c(:,3)-9856.*c(:,1).^2.*c(:,2).^3.*c(:,4)+1344.*c(:,1).^2.*c(:,2).^2.*(-11.*c(:,3).^2+8.*c(:,1).*c(:,5))-3072.*c(:,1).^3.*c(:,2).*(-7.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))+512.*c(:,1).^3.*(7.*c(:,3).^3-24.*c(:,1).*c(:,3).*c(:,5)+4.*c(:,1).*(-3.*c(:,4).^2+8.*c(:,1).*c(:,7)))), ...
    c1inv.^(27/4)/(262144).*(4807.*c(:,2).^7-35112.*c(:,1).*c(:,2).^5.*c(:,3)+36960.*c(:,1).^2.*c(:,2).^4.*c(:,4)-4928.*c(:,1).^2.*c(:,2).^3.*(-15.*c(:,3).^2+8.*c(:,1).*c(:,5))+10752.*c(:,1).^3.*c(:,2).^2.*(-11.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))-512.*c(:,1).^3.*c(:,2).*(77.*c(:,3).^3-168.*c(:,1).*c(:,3).*c(:,5)+12.*c(:,1).*(-7.*c(:,4).^2+8.*c(:,1).*c(:,7)))+2048.*c(:,1).^4.*(21.*c(:,3).^2.*c(:,4)-24.*c(:,1).*c(:,3).*c(:,6)+8.*c(:,1).*(-3.*c(:,4).*c(:,5)+4.*c(:,1).*c(:,8)))), ...
    c1inv.^(31/4)/(8388608).*(-129789.*c(:,2).^8+1076768.*c(:,1).*c(:,2).^6.*c(:,3)-1123584.*c(:,1).^2.*c(:,2).^5.*c(:,4)+147840.*c(:,1).^2.*c(:,2).^4.*(-19.*c(:,3).^2+8.*c(:,1).*c(:,5))-315392.*c(:,1).^3.*c(:,2).^3.*(-15.*c(:,3).*c(:,4)+4.*c(:,1).*c(:,6))+43008.*c(:,1).^3.*c(:,2).^2.*(55.*c(:,3).^3-88.*c(:,1).*c(:,3).*c(:,5)+4.*c(:,1).*(-11.*c(:,4).^2+8.*c(:,1).*c(:,7)))-49152.*c(:,1).^4.*c(:,2).*(77.*c(:,3).^2.*c(:,4)-56.*c(:,1).*c(:,3).*c(:,6)+8.*c(:,1).*(-7.*c(:,4).*c(:,5)+4.*c(:,1).*c(:,8)))+4096.*c(:,1).^4.*(-77.*c(:,3).^4+336.*c(:,1).*c(:,3).^2.*c(:,5)-48.*c(:,1).*c(:,3).*(-7.*c(:,4).^2+8.*c(:,1).*c(:,7))+64.*c(:,1).^2.*(-3.*c(:,5).^2-6.*c(:,4).*c(:,6)+8.*c(:,1).*c(:,9))))];


% A = [1./a(:,1),...
%     -(a(:,2)./a(:,1).^3),...
%     (2.*a(:,2).^2-a(:,1).*a(:,3))./a(:,1).^5,...
%     -((5.*a(:,2).^3-5.*a(:,1).*a(:,2).*a(:,3)+a(:,1).^2.*a(:,4))./a(:,1).^7),...
%     1./a(:,1).^9.*(14.*a(:,2).^4-21.*a(:,1).*a(:,2).^2.*a(:,3)+3.*a(:,1).^2.*a(:,3).^2+6.*a(:,1).^2.*a(:,2).*a(:,4)-a(:,1).^3.*a(:,5)),...
%     1./a(:,1).^11.*(-42.*a(:,2).^5+84.*a(:,1).*a(:,2).^3.*a(:,3)-28.*a(:,1).^2.*a(:,2).^2.*a(:,4)+7.*a(:,1).^2.*a(:,2).*(-4.*a(:,3).^2+a(:,1).*a(:,5))+a(:,1).^3.*(7.*a(:,3).*a(:,4)-a(:,1).*a(:,6))),...
%     1./a(:,1).^13.*(132.*a(:,2).^6-330.*a(:,1).*a(:,2).^4.*a(:,3)+120.*a(:,1).^2.*a(:,2).^3.*a(:,4)-36.*a(:,1).^2.*a(:,2).^2.*(-5.*a(:,3).^2+a(:,1).*a(:,5))+8.*a(:,1).^3.*a(:,2).*(-9.*a(:,3).*a(:,4)+a(:,1).*a(:,6))+a(:,1).^3.*(-12.*a(:,3).^3+8.*a(:,1).*a(:,3).*a(:,5)+a(:,1).*(4.*a(:,4).^2-a(:,1).*a(:,7)))),...
%     1./a(:,1).^15.*(-429.*a(:,2).^7+1287.*a(:,1).*a(:,2).^5.*a(:,3)-495.*a(:,1).^2.*a(:,2).^4.*a(:,4)+165.*a(:,1).^2.*a(:,2).^3.*(-6.*a(:,3).^2+a(:,1).*a(:,5))-45.*a(:,1).^3.*a(:,2).^2.*(-11.*a(:,3).*a(:,4)+a(:,1).*a(:,6))+3.*a(:,1).^3.*a(:,2).*(55.*a(:,3).^3-30.*a(:,1).*a(:,3).*a(:,5)+3.*a(:,1).*(-5.*a(:,4).^2+a(:,1).*a(:,7)))+a(:,1).^4.*(-45.*a(:,3).^2.*a(:,4)+9.*a(:,1).*a(:,3).*a(:,6)+a(:,1).*(9.*a(:,4).*a(:,5)-a(:,1).*a(:,8)))),...
%     1/a(:,1).^17.*(1430.*a(:,2).^8-5005.*a(:,1).*a(:,2).^6.*a(:,3)+2002.*a(:,1).^2.*a(:,2).^5.*a(:,4)-715.*a(:,1).^2.*a(:,2).^4.*(-7.*a(:,3).^2+a(:,1).*a(:,5))+220.*a(:,1).^3.*a(:,2).^3.*(-13.*a(:,3).*a(:,4)+a(:,1).*a(:,6))-55.*a(:,1).^3.*a(:,2).^2.*(26.*a(:,3).^3-12.*a(:,1).*a(:,3).*a(:,5)+a(:,1).*(-6.*a(:,4).^2+a(:,1).*a(:,7)))+10.*a(:,1).^4.*a(:,2).*(66.*a(:,3).^2.*a(:,4)-11.*a(:,1).*a(:,3).*a(:,6)+a(:,1).*(-11.*a(:,4).*a(:,5)+a(:,1).*a(:,8)))+a(:,1).^4.*(55.*a(:,3).^4-55.*a(:,1).*a(:,3).^2.*a(:,5)+5.*a(:,1).*a(:,3).*(-11.*a(:,4).^2+2.*a(:,1).*a(:,7))+a(:,1).^2.*(5.*a(:,5).^2+10.*a(:,4).*a(:,6)-a(:,1).*a(:,9))))];

A = [1.*c1inv.^(1/4),...
    -(a(:,2).*c1inv.^(3/4)),...
    (2.*a(:,2).^2-a(:,1).*a(:,3)).*c1inv.^(5/4),...
    -((5.*a(:,2).^3-5.*a(:,1).*a(:,2).*a(:,3)+a(:,1).^2.*a(:,4)).*c1inv.^(7/4)),...
    1.*c1inv.^(9/4).*(14.*a(:,2).^4-21.*a(:,1).*a(:,2).^2.*a(:,3)+3.*a(:,1).^2.*a(:,3).^2+6.*a(:,1).^2.*a(:,2).*a(:,4)-a(:,1).^3.*a(:,5)),...
    1.*c1inv.^(11/4).*(-42.*a(:,2).^5+84.*a(:,1).*a(:,2).^3.*a(:,3)-28.*a(:,1).^2.*a(:,2).^2.*a(:,4)+7.*a(:,1).^2.*a(:,2).*(-4.*a(:,3).^2+a(:,1).*a(:,5))+a(:,1).^3.*(7.*a(:,3).*a(:,4)-a(:,1).*a(:,6))),...
    1.*c1inv.^(13/4).*(132.*a(:,2).^6-330.*a(:,1).*a(:,2).^4.*a(:,3)+120.*a(:,1).^2.*a(:,2).^3.*a(:,4)-36.*a(:,1).^2.*a(:,2).^2.*(-5.*a(:,3).^2+a(:,1).*a(:,5))+8.*a(:,1).^3.*a(:,2).*(-9.*a(:,3).*a(:,4)+a(:,1).*a(:,6))+a(:,1).^3.*(-12.*a(:,3).^3+8.*a(:,1).*a(:,3).*a(:,5)+a(:,1).*(4.*a(:,4).^2-a(:,1).*a(:,7)))),...
    1.*c1inv.^(15/4).*(-429.*a(:,2).^7+1287.*a(:,1).*a(:,2).^5.*a(:,3)-495.*a(:,1).^2.*a(:,2).^4.*a(:,4)+165.*a(:,1).^2.*a(:,2).^3.*(-6.*a(:,3).^2+a(:,1).*a(:,5))-45.*a(:,1).^3.*a(:,2).^2.*(-11.*a(:,3).*a(:,4)+a(:,1).*a(:,6))+3.*a(:,1).^3.*a(:,2).*(55.*a(:,3).^3-30.*a(:,1).*a(:,3).*a(:,5)+3.*a(:,1).*(-5.*a(:,4).^2+a(:,1).*a(:,7)))+a(:,1).^4.*(-45.*a(:,3).^2.*a(:,4)+9.*a(:,1).*a(:,3).*a(:,6)+a(:,1).*(9.*a(:,4).*a(:,5)-a(:,1).*a(:,8)))),...
    1*c1inv.^(17/4).*(1430.*a(:,2).^8-5005.*a(:,1).*a(:,2).^6.*a(:,3)+2002.*a(:,1).^2.*a(:,2).^5.*a(:,4)-715.*a(:,1).^2.*a(:,2).^4.*(-7.*a(:,3).^2+a(:,1).*a(:,5))+220.*a(:,1).^3.*a(:,2).^3.*(-13.*a(:,3).*a(:,4)+a(:,1).*a(:,6))-55.*a(:,1).^3.*a(:,2).^2.*(26.*a(:,3).^3-12.*a(:,1).*a(:,3).*a(:,5)+a(:,1).*(-6.*a(:,4).^2+a(:,1).*a(:,7)))+10.*a(:,1).^4.*a(:,2).*(66.*a(:,3).^2.*a(:,4)-11.*a(:,1).*a(:,3).*a(:,6)+a(:,1).*(-11.*a(:,4).*a(:,5)+a(:,1).*a(:,8)))+a(:,1).^4.*(55.*a(:,3).^4-55.*a(:,1).*a(:,3).^2.*a(:,5)+5.*a(:,1).*a(:,3).*(-11.*a(:,4).^2+2.*a(:,1).*a(:,7))+a(:,1).^2.*(5.*a(:,5).^2+10.*a(:,4).*a(:,6)-a(:,1).*a(:,9))))];

Tscal = cell(1,nvp.degree);
TscalInv = cell(1,nvp.degree);
TscalInv{1} = invertibleMatrix(diag(a(:,1)),diag(A(:,1)));
Tscal{1} = invertibleMatrix(diag(A(:,1)),diag(a(:,1)));
for i=2:nvp.degree
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

function [sigmaSquared, TinOd, vbar, wbar] = inputNormalOutputDiagonalTransformation(v, w, nvp)
%inputNormalOutputDiagonalTransformation Return a polynomial input-normal/output-diagonal transformation x = Φ(z).
%
%   Usage: [sigmaSquared,Tbar] = inputNormalOutputDiagonalTransformation(v, w)
%
%   Inputs:
%       v,w     - cell arrays containing the polynomial energy function coefficients.
%
%     Optional name/value pair inputs:
%        degree - desired degree of the balanced realization. Ex. for linear
%                 dynamics, choosing degree=1 will produce quadratic
%                 approximations for the energy functions, linear approximations
%                 for the transformations, and the output will be a linear
%                 balanced realization. The default will be on less than the
%                 degree of the v, i.e. we will balance everything available, no
%                 more and no less.
%       r       - reduced-order dimension if computing a balance-and-reduce
%                 transformation.
%       verbose - optional argument to print runtime information
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
    nvp.degree = length(v) - 1
    nvp.r = sqrt(numel(v{2}))
    nvp.verbose = false
end
vec = @(X) X(:); % Create a vec function for readability
if nvp.verbose
    fprintf('Computing the degree %d input-normal/output-diagonal balancing transformation...\n', nvp.degree)
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

% Approach 1: Full-order transformation 
% Tin = invertibleMatrix(Rinv.'*V,  Xi\U.'*L.');

% Approach 2: Reduced-order transformation
% Xir = Xi(1:nvp.r,1:nvp.r); Ur = U(:,1:nvp.r); Vr = V(:,1:nvp.r);
% Tinr = invertibleMatrix(Rinv.'*Vr,  Xir\Ur.'*L.');

% Approach 3: Full sized but reduced-order transformation
XirInv = diag([diag(1./Xi(1:nvp.r,1:nvp.r));zeros(n-nvp.r,1)]);
U = [U(:,1:nvp.r), zeros(n,n-nvp.r)]; V = [V(:,1:nvp.r), zeros(n,n-nvp.r)];
Tin = invertibleMatrix(Rinv.'*V,  XirInv*U.'*L.');



%%%
% Possibly use a hybrid of Tin with the full forward but reduced inverse
% transformation
%%%

% Transform the energy functions using Tin; 
% ** this does not involve any inversions of Tin, so we can use full Tin **
[vtilde, wtilde] = transformEnergyFunctionsLinear(v, w, Tin);

% Name Ṽ₂ and W̃₂; in principle they would be the first two entries, but
% analytically we know what they are
% V2tilde = reshape(vtilde{2},n,n); W2tilde = reshape(wtilde{2},n,n);
V2tilde = speye(n); vtilde{2} = vec(V2tilde);
W2tilde = sparse(diag(diag(Xi))) .^ 2; wtilde{2} = vec(W2tilde);

%% Step 2: Compute the higher-order terms in the second transformation
% Preallocate the cell array, the first term is identity
Tod = cell(1, nvp.degree);
Tod{1} = speye(n);

retainedStates = (1:n) <= nvp.r;
entries2keep = retainedStates;             
% entries2keep = vec(vec(entries2keep) & retainedStates);

% Compute the higher-order terms according to Corollary 1 [1]
for k = 3:nvp.degree + 1
    if nvp.verbose; fprintf("    Computing degree %i coefficient... ", k - 1); tic; end
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
        if i > length(vtilde); break; end
        temp = temp + calTTv(Tod, i, k, vtilde{i}); % TODO: Accelerate this
    end
    RHS = [RHS; -Nk * temp];
    
    %% Form output-diagonal equations coefficient matrix
    CoeffMatrix = [CoeffMatrix; 2 * Nkhat * kron(speye(n^(k-1)), W2tilde)]; % TODO: kronecker rules; currently this is a sparse repmat diag kind of thing, so it is ok
    
    temp = zeros(size(Nkhat, 2), 1);
    for i = 2:k - 2
        j = k - i;
        temp = temp + vec(Tod{j}.' * W2tilde * Tod{i});
    end
    for i = 3:k
        if i > length(wtilde); break; end
        temp = temp + calTTv(Tod, i, k, wtilde{i}); % TODO: Accelerate this; should already be efficient via calTTv
    end
    RHS = [RHS; -Nkhat * temp];
    
    %% Form `flexibility' equations (Kronecker product repeated entries)
    [linclassidx] = referenceElementMap(n, k-1);
    
    linclassidx(linclassidx) = []; % Basically remove the reference element so one is nonzero and the rest we eliminate
    idxs = vec(( n*(linclassidx-1) + (1:n) ).');
    
    %% Set extra parameter equation
    % TODO: find best parameters; for now just solve "a" solution with mldivide
    % parameterEqsRHS = 1;
    % parameterEqsCoeff = zeros(1,n^k);
    % parameterEqsCoeff(4) = 1;
    
    % CoeffMatrix(:,4) = [];
    % indices(4) = [];
    
    % idxs = [idxs; 16]; % alternatively could solve minimum norm solution
    
    method = 'full-order';
    switch method
        case 'full-order'
            %% Method 1: Solve equations (full-order)
            % Assemble equations
            CoeffMatrix(:, idxs) = []; % Kronecker flexibility

            % Form index set for the nonzero transformation components
            indices = 1:n^k; indices(idxs) = [];

            Tod{k-1} = zeros(n, n^(k-1));
            Tod{k-1}(indices) = CoeffMatrix \ RHS;                     % Method 1: matlab uses sparse QR from SuiteSparseQR
            % Tod{k - 1}(indices) = lsqminnorm(CoeffMatrix, RHS);      % Method 2: minimum norm solution

            % Now "truncate" -> set to zero
            % entries2keep = vec(vec(entries2keep) & retainedStates);
            % Tod{k-1}(:, ~entries2keep) = 0;

        case 'reduced-order'
            %% Method 2: Solve equations (reduced-order)
            % Assemble equations
            CoeffMatrix(:, idxs) = []; % Kronecker flexibility

            % Form index set for the nonzero transformation components
            indices = 1:n^k; indices(idxs) = [];

            Tod{k-1} = zeros(n, n^(k-1));
            Tod{k-1}(indices) = CoeffMatrix \ RHS;                     % Method 1: matlab uses sparse QR from SuiteSparseQR
            % Tod{k - 1}(indices) = lsqminnorm(CoeffMatrix, RHS);      % Method 2: minimum norm solution

            entries2keep = vec(vec(entries2keep) & retainedStates);
            Tod{k-1} = Tod{k-1}(:, entries2keep);

    end
    if nvp.verbose; fprintf("completed in %f seconds. \n", toc); end
end

%% Combine transformation with linear input-normal transformation
TinOd = cell(1,nvp.degree);
TinOd{1} = Tin;
for k = 2:nvp.degree
    TinOd{k} = Tin * Tod{k};
    TinOd{k} = kronMonomialSymmetrize(TinOd{k}, n, k); % Symmetrize the transformation rows
end

%% Pluck out the singular value function coefficients
[vbar, wbar] = transformEnergyFunctions(v, w, TinOd, true); % Could transform just the observability; could probably even just compute the diagonal entries

sigmaSquared = zeros(n, nvp.degree);
for k = 2:nvp.degree + 1
    if k > length(wbar); break; end
    if nvp.verbose
        [N] = equivalenceClassIndices(n, k);
        
        fprintf("      - The largest entry in v%i is %.1e; ", k, max(abs(N * vbar{k}))) % Should be zero, other than the first time which is one
        fprintf("the largest off-diagonal entry in w%i is %.1e\n", k, max(abs(N(n + 1:end, :) * wbar{k}))) % Should be diagonal
        
        sigmaSquared(:, k - 1) = N(1:n, :) * wbar{k}; % Since the index set is already computed
    else
        indexSet = linspace(1, n^k, n);
        sigmaSquared(:, k-1) = wbar{k}(indexSet);
    end
end

if nvp.verbose
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
