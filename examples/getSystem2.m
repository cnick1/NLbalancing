function [f, g, h] = getSystem2(kawanoModel)
%getSystem2  Generates a 2D polynomial system for testing energy functions.
%            Model based on systems taken from [1,2].
%
%   Usage:  [f,g,h] = getSystem2()
%
%   Inputs:
%         kawanoModel - whether to use the model from [2] (true)
%                ẋ₁ = -x₁ + x₂ - x₂² + u₁ + 2 x₂ u₁ - 0.05 x₁ x₂ u₁
%                ẋ₂ =     - x₂       + u₁           - 0.05 x₂² u₁
%                 y =  x₁
%         or the approximation used in [1]
%                ẋ₁ = -x₁ + x₂ - x₂² + u₁ + 2 x₂ u₁ - 0.05 x₁ x₂ u₁
%                ẋ₂ =     - x₂       + u₁           - 0.05 x₂² u₁
%                 y =  x₁ + x₂
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Note: In [1], the higher-order g terms are neglected. [3] uses
%   the true model from [2].
%
%   Reference: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki,
%               “Scalable computation of energy functions for nonlinear
%               balanced truncation,” Computer Methods in Applied Mechanics
%               and Engineering, vol. 427, p. 117011, Jul. 2024, doi:
%               10.1016/j.cma.2024.117011
%              [2] Y. Kawano and J. M. A. Scherpen, “Model reduction by
%               differential balancing based on nonlinear hankel operators,”
%               IEEE Transactions on Automatic Control, vol. 62, no. 7,
%               pp. 3293–3308, Jul. 2017, doi: 10.1109/tac.2016.2628201
%              [3] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗∞
%               energy functions for polynomial control-affine systems,"
%               IEEE Transactions on Automatic Control, pp. 1–13, 2024,
%               doi: 10.1109/tac.2024.3494472
%%

if nargin < 1
     kawanoModel = false;
end

if kawanoModel % Use Kawano model
     A = [-1 1;
          0 -1];
     N = [0 0 0 -1;
          0 0 0 0];
     B = [1;
          1];
     C = [1 0]; % Main difference with Kawano model
     G1 = [0 2;
          0 0];
     %   G2 = [0 -0.05 0 0; % G2 is not in Kawano model, so ignore if desired
     %         0 0 0 -0.05];
else % Use modified model from Kramer et. al. (default)
     A = [-1 1;
          0 -1];
     N = [0 0 0 -1;
          0 0 0 0];
     B = [1;
          1];
     C = [1 1];
     G1 = [0 2;
          0 0];
     %   G2 = [0 -0.05 0 0;
     %         0 0 0 -0.05];
end

f = {A, N};
g = {B, G1}; %, G2};
h = {C};

end
