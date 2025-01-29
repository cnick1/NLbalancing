function [f, g, h] = getSystem7()
%getSystem7  Generates an aircraft stall model from Garrard 1977 [1]
%  for testing stabilization.
%
%   Usage:  [f,g,h] = getSystem7()
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The dynamics correspond to the input-output system from [1]
%         ẋ₁ = x₃ - x₁² x₃ - 0.088 x₁ x₃ - 0.877 x₁ + 0.47 x₁² + 3.846 x₁³
%              - 0.215 u + 0.28 u x₁² + 0.47 u² x₁ +  0.63 u³ - 0.019 x₂²
%         ẋ₂ = x₃
%         ẋ₃ = -0.396 x₃ - 4.208 x₁ - 0.47 x₁² - 3.564 x₁³ - 20.967 u +
%              6.265 ux₁² + 46 u² + 61.4 u³
%
%   Reference: [1] W. L. Garrard and J. M. Jordan, “Design of nonlinear
%               automatic flight control systems,” Automatica, vol. 13,
%               no. 5, pp. 497–505, Sep. 1977,
%               doi: 10.1016/0005-1098(77)90070-x
%
%   Part of the NLbalancing repository.
%%

A = [-0.877 0 1;
     0 0 1;
     -4.208 0 -0.396];
B = [-0.215;
     0;
     -20.967];
C = .5 * eye(3);
F2 = sparse(3, 9); F2(1, 1) = 0.47; F2(1, 3) = -0.088; F2(1, 5) = -0.019; F2(3, 1) = -0.47;
F3 = sparse(3, 27); F3(1, 1) = 3.846; F3(1, 3) = -1; F3(3, 1) = -3.564;
G1 = sparse(3, 3);
G2 = sparse(3, 9); G2(1, 1) = 0.28; G2(3, 1) = 6.265;

f = {A, F2, F3};
g = {B, G1, G2};
h = {C};

end
