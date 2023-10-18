function [f, g, h] = getSystem7()
%getSystem7  Generates an aircraft stall model from Garrard 1977 [1]
%  for testing stabilization.
%
%   Usage:  [f,g,h] = getSystem7()
%
%   The "matrices" correspond to the input-output system
%
%         \dot{x}_1 & = x_3 - x_1^2 x_3 - 0.088 x_1 x_3 - 0.877 x_1 + 0.47 x_1^2 + 3.846 x_1^3 - 0.215 u + 0.28 u x_1^2 + 0.47 u^2 x_1 +  0.63 u^3 - 0.019 x_2^2 \\
%         \dot{x}_2 & = x_3                                                                                                                                      \\
%         \dot{x}_3 & = -0.396 x_3 - 4.208 x_1 - 0.47 x_1^2 - 3.564 x_1^3 - 20.967 u + 6.265 ux_1^2 + 46 u^2 + 61.4 u^3
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
