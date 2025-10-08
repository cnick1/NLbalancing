function [f, g, h] = getSystem25(lienard, eps)
%getSystem25 Generates a stable Van der Pol oscillator for testing
%            computing the region of attraction using Zubov's method.
%
%   Usage:  [f,g,h] = getSystem25()
%
%   Inputs:    lienard - whether to use model (1) or (2)
%              eps     - damping parameter
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The Van der Pol oscillator is a nonlinear oscillator with
%   nonlinear damping. The basic equation in second-order form is
%
%               ẍ + ϵ(1 - x²)ẋ + x = u(t)
%
%   This can be put in first-order form with x₁ = x, x₂ = ẋ as   
%
%             ẋ₁ = x₂
%      (1)    ẋ₂ = -x₁ - ϵ(1 - x₁²)x₂ + u(t)
%              y = x₁
% 
%   Alternatively, [3] defines z₁ = x, z₂ = ẋ + ϵ(x - x³/3) to obtain
% 
%             ̇z₁ = z₂ - ϵ(z₁ - z₁³/3)
%      (2)    ̇z₂ = -z₁
%              y = z₁
% 
%   Is it possible to form a coordinate transformation (bijection) between
%   (1) and (2) without going to second-order? YES! The transformation is
%
%               x = Φ(z)   = [z₁; z₂ - ϵ(z₁ - z₁³/3)]
%               z = Φ⁻¹(x) = [x₁; x₂ + ϵ(x₁ - x₁³/3)]
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   !!!!!!!!!! This is sort of like a near identity transformation !!!!!!!!
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
%   This is similar to Lienard transformation x₂ = x - x³/3 + ẋ/ϵ, with
%   which the system can be put in the form
%
%             ẋ₁ = -ϵ(x₁ - x₂ - x₁³/3)
%      (3)    ẋ₂ = -x₁/ϵ
%              y = x₁
%
%   Things I am curious about: what do the Lyapunov functions, regions of
%   attraction, etc. look like in these different coordinate systems? Is
%   one more suitable for polynomial approximations? How does this
%   transformation relate to a nonlinear balancing transformation?
%
%   References: [1] https://en.wikipedia.org/wiki/Van_der_Pol_oscillator
%               [2] H. K. Khalil, Nonlinear systems, Third edition, Pearson
%                Education, 2013
%               [3] S. Margolis and W. Vogt, "Control engineering
%                applications of V. I. Zubov's construction procedure for
%                Lyapunov functions," IEEE Transactions on Automatic
%                Control, vol. 8, no. 2, pp. 104–113, Apr. 1963, doi:
%                10.1109/tac.1963.1105553
%
%
%   Notes:
%    - The standard form in [1] is for an unstable Van der Pol
%      oscillator, which corresponds to μ = -ϵ.
%    - In [2], instead of considering the stable Van der Pol oscillator,
%      Khalil considers the unstable oscillator integrated backward in
%      time. The result is similar, but I prefer the physical intuition of
%      the negative damping.
%
%%
if nargin < 2
    eps = 1;
    if nargin < 1
        lienard = false;
    end
end

n = 2;

if lienard % form from [3]
    %       ẋ₁ = x₂ - ϵ(x₁ - x₁³/3)
    %  (2)  ẋ₂ = -x₁
    %        y =  x₁

    A = [-eps 1;
        -1 0];
    F2 = zeros(n,n^2);
    F3 = zeros(n,n^3); F3(1,1) = eps/3;

    B = [0; 1];

else % standard state-space form
    %       ẋ₁ = x₂
    %  (1)  ẋ₂ = -x₁ - ϵ(1 - x₁²)x₂ + u(t)
    %        y =  x₁

    A = [0 1;
        -1 -eps];
    F2 = zeros(n,n^2);
    F3 = zeros(n,n^3); F3(2,2) = eps;

    B = [0; 1];

end

C = [1; 0];

f = {A, F2, F3};
g = {B};
h = {C};

end
