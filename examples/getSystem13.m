function [f, g, h] = getSystem13()
%getSystem13  Generates a polynomial approximation to the 2D model from Gray and Scherpen 2001 [1]
%
%   Usage:  [f,g,h] = getSystem13()
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The output is a polynomial approximation to the system given by
%           f(x) = -[α² x₁ + 2 α x₂ + (α² - 2)x₂²;
%                                x₂              ]
%           g(x) = √2[α - 2 x₂;
%                        1    ]
%           h(x) = 1/√3 [3 α (x₁ + x₂²) + (α - 2√2)x₂]
%
%    where α = (√3 + √2)(√3 + 2)
%
%   References: [1] W. S. Gray and J. M. A. Scherpen, “On the nonuniqueness
%               of singular value functions and balanced nonlinear
%               realizations,” Systems & Control Letters, vol. 44, no. 3,
%               pp. 219–232, Oct. 2001, doi: 10.1016/s0167-6911(01)00144-x
%
%%

degree = 2;

n = 2; x = sym('x', [1, n]).'; syms(x);

alpha = (sqrt(3) + sqrt(2)) * (sqrt(3) + 2);

fsym = - [alpha ^ 2 * x1 + 2 * alpha * x2 + (alpha ^ 2 - 2) * x2 ^ 2;
    x2];
gsym = sqrt(2) * [alpha - 2 * x2;
    1];
hsym = 1 / sqrt(3) * (3 * alpha * (x1 + x2 ^ 2) + (alpha - 2 * sqrt(2)) * x2);

[f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, degree);

f{1} = full(f{1}); g{1} = full(g{1}); h{1} = full(h{1});
f{2} = full(f{2}); g{2} = full(g{2}); h{2} = full(h{2});

end
