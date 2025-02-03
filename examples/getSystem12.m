function [f, g, h] = getSystem12(degree, transformedModel)
%getSystem12  Generates a polynomial approximation to the 2D model from Fujimoto and Scherpen 2001, 2005, 2010 [1-3]
%
%   Usage:  [f,g,h] = getSystem12(degree, transformedModel)
%
%   Inputs:
%                  degree  -  degree d taylor approximation
%        transformedModel  -  whether or not to use the transformed (input-normal) model
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The model from [1-3] is
%       f(x) = [-9x₁ + 6x₁² x₂ + 6x₂³ - x₁⁵ - 2x₁³ x₂² - x₁ x₂⁴
%               -9x₂ - 6x₁³ - 6x₁x₂² - x₁⁴ x₂ - 2x₁² x₂³ - x₂⁵]
%       g(x) = [\frac{3√2(9-6x₁x₂+x₁⁴-x₂⁴)}{9+x₁⁴+2x₁²x₂²+x₂⁴},
%                 \frac{√2(-9x₁² - 27 x₂² + 6 x₁³ x₂ + 6 x₁ x₂³ - (x₁² + x₂²)³)}{9+x₁⁴+2x₁²x₂²+x₂⁴};
%               \frac{√2(27x₁²+9x₂²+6x₁³x₂+6x₁x₂³+(x₁²+x₂²)³}{9+x₁⁴+2x₁²x₂²+x₂⁴},
%                 \frac{3√2(9 + 6 x₁ x₂  - x₁⁴ + x₂⁴)}{9+x₁⁴+2x₁²x₂²+x₂⁴}]
%       h(x) = [\frac{2√2(3x₁ + x₁ x₂² + x₂³)(3 - x₁⁴ - 2x₁² x₂² - x₂⁴)}{1 + x₁⁴ + 2 x₁² x₂² + x₂⁴};
%                 \frac{√2(3x₂ - x₁³ - x₁ x₂²)(3 - x₁⁴ - 2 x₁² x₂² - x₂⁴)}{1 + x₁⁴ + 2 x₁² x₂² + x₂⁴}]
%
%   References: [1] K. Fujimoto and J. M. A. Scherpen, “Model reduction
%                for nonlinear systems based on the differential
%                eigenstructure of Hankel operators,” in Proceedings of
%                the 40th IEEE Conference on Decision and Control (Cat.
%                No.01CH37228), IEEE, 2001. doi: 10.1109/cdc.2001.980322
%               [2] K. Fujimoto and J. M. A. Scherpen, “Nonlinear
%                input-normal realizations based on the differential
%                eigenstructure of Hankel operators,” IEEE Transactions
%                on Automatic Control, vol. 50, no. 1, pp. 2–18, Jan.
%                2005, doi: 10.1109/tac.2004.840476
%               [3] K. Fujimoto and J. M. A. Scherpen, “Balanced
%                realization and model order reduction for nonlinear
%                systems based on singular value analysis,” SIAM Journal
%                on Control and Optimization, vol. 48, no. 7, pp.
%                4591–4623, Jan. 2010, doi: 10.1137/070695332
%
%%

if nargin < 2
        transformedModel = true;
        if nargin < 1
                degree = 9;
        end
end

n = 2;
x = sym('x', [1, n]).';
syms(x);

if transformedModel % Use the model in the transformed coordinates (21) in Example 3 Fujimoto 2001
        fsym = [-9 * x1 - x1 ^ 5 - 2 * x1 ^ 3 * x2 ^ 2 - x1 * x2 ^ 4;
                -9 * x2 - x1 ^ 4 * x2 - 2 * x1 ^ 2 * x2 ^ 3 - x2 ^ 5];
        gsym = [sqrt(18 + 2 * x1 ^ 4 + 4 * x1 ^ 2 * x2 ^ 2 + 2 * x2 ^ 4), 0;
                0, sqrt(18 + 2 * x1 ^ 4 + 4 * x1 ^ 2 * x2 ^ 2 + 2 * x2 ^ 4)];
        hsym = [(6 * x1 - 2 * x1 ^ 5 - 4 * x1 ^ 3 * x2 ^ 2 - 2 * x1 * x2 ^ 4) * sqrt(18 + 2 * x1 ^ 4 + 4 * x1 ^ 2 * x2 ^ 2 + 2 * x2 ^ 4) / (1 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4);
                (3 * x2 - x1 ^ 4 * x2 - 2 * x1 ^ 2 * x2 ^ 3 -  x2 ^ 5) * sqrt(18 + 2 * x1 ^ 4 + 4 * x1 ^ 2 * x2 ^ 2 + 2 * x2 ^ 4) / (1 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4)];
else % Use the model in the original coordinates (Example 1) in Fujimoto 2001
        fsym = [-9 * x1 + 6 * x1 ^ 2 * x2 + 6 * x2 ^ 3 - x1 ^ 5 - 2 * x1 ^ 3 * x2 ^ 2 - x1 * x2 ^ 4;
                -9 * x2 - 6 * x1 ^ 3 - 6 * x1 * x2 ^ 2 - x1 ^ 4 * x2 - 2 * x1 ^ 2 * x2 ^ 3 - x2 ^ 5];
        gsym = [3 * sqrt(2) * (9 - 6 * x1 * x2 + x1 ^ 4 - x2 ^ 4) / (9 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4), sqrt(2) * (-9 * x1 ^ 2 - 27 * x2 ^ 2 + 6 * x1 ^ 3 * x2 + 6 * x1 * x2 ^ 3 - (x1 ^ 2 + x2 ^ 2) ^ 3) / (9 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4);
                sqrt(2) * (27 * x1 ^ 2 + 9 * x2 ^ 2 + 6 * x1 ^ 3 * x2 + 6 * x1 * x2 ^ 3 + (x1 ^ 2 + x2 ^ 2) ^ 3) / (9 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4), 3 * sqrt(2) * (9 + 6 * x1 * x2 - x1 ^ 4 + x2 ^ 4) / (9 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4)];
        % hsym = [2 * sqrt(2) * (3 * x1 + x1 * x2 ^ 2 + x2 ^ 3) * (3 - x1 ^ 4 - 2 * x1 ^ 2 * x2 ^ 2 - x2 ^ 4) / (1 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4); % This is what's in 2001 paper, but it has a typo and so doesn't match the results they present
        hsym = [2 * sqrt(2) * (3 * x1 + x1 ^ 2 * x2 + x2 ^ 3) * (3 - x1 ^ 4 - 2 * x1 ^ 2 * x2 ^ 2 - x2 ^ 4) / (1 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4);
                sqrt(2) * (3 * x2 - x1 ^ 3 - x1 * x2 ^ 2) * (3 - x1 ^ 4 - 2 * x1 ^ 2 * x2 ^ 2 - x2 ^ 4) / (1 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4)];
end
[f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, degree);

end
