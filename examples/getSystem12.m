function [f, g, h] = getSystem12(degree, transformedModel)
%getSystem12  Generates a polynomial approximation to the 2D model from Fujimoto and Scherpen 2001, 2005, 2010 [1-3]
%
%   Usage:  [f,g,h] = getSystem12(degree)
%
%       f(x) & = \begin{pmatrix}
%       -9x_1 + 6x_1^2 x_2 + 6x_2^3 - x_1^5 - 2x_1^3 x_2^2 - x_1 x_2^4 \\
%       -9x_2 - 6x_1^3 - 6x_1x_2^2 - x_1^4 x_2 - 2x_1^2 x_2^3 - x_2^5
%       \end{pmatrix} \\
%       g(x) & = \begin{pmatrix}
%           \frac{3\sqrt{2}(9-6x_1x_2+x_1^4-x_2^4)}{9+x_1^4+2x_1^2x_2^2+x_2^4} & \frac{\sqrt{2}(-9x_1^2 - 27 x_2^2 + 6 x_1^3 x_2 + 6 x_1 x_2^3 - (x_1^2 + x_2^2)^3)}{9+x_1^4+2x_1^2x_2^2+x_2^4} \\
%           \frac{\sqrt{2}(27x_1^2+9x_2^2+6x_1^3x_2+6x_1x_2^3+(x_1^2+x_2^2)^3}{9+x_1^4+2x_1^2x_2^2+x_2^4} & \frac{3\sqrt{2}(9 + 6 x_1 x_2  - x_1^4 + x_2^4)}{9+x_1^4+2x_1^2x_2^2+x_2^4}
%       \end{pmatrix} \\
%       h(x) & = \begin{pmatrix}
%       \frac{2\sqrt{2}(3x_1 + x_1 x_2^2 + x_2^3)(3 - x_1^4 - 2x_1^2 x_2^2 - x_2^4)}{1 + x_1^4 + 2 x_1^2 x_2^2 + x_2^4} \\
%       \frac{\sqrt{2}(3x_2 - x_1^3 - x_1 x_2^2)(3 - x_1^4 - 2 x_1^2 x_2^2 - x_2^4)}{1 + x_1^4 + 2 x_1^2 x_2^2 + x_2^4}
%       \end{pmatrix}
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
else % Use the model in the transformed coordinates (Example 1) in Fujimoto 2001
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
