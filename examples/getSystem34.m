function [f, g, h] = getSystem34()
%getSystem34 Returns Hirsch's 2D bio model
%
%   Usage:  [f,g,h] = getSystem34()
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The dynamics correspond to the input-output system
%
%           ẋ₁ = −x₁ + λ x₁x₂ + u₁,
%           ẋ₂ = −x₂ + λ x₁x₂ + u₂,
%            y = [x₁,x₂].
%
%   References: [1] 
%
%%
n = 3; x = sym('x', [1, n]).';
lambda = 0.1;
fsym = [-x(1) + lambda * x(1)*x(2);
    -0.1* x(2) + lambda * x(1)*x(3);
    -0.1* x(3) + x(1) + lambda * x(2)*x(3);
];
gsym = [1 1; 0 0.1; 0 0];
hsym = x;

[f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, 3);

end
