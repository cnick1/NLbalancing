function [f, g, h] = getSystem1()
%getSystem1  Generates a simple 1D polynomial system for testing energy functions.
%
%   Usage:  [f, g, h] = getSystem1()
%
%   Outputs:    f,g,h  - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The dynamics correspond to the input-output system
%
%       ẋ = -2x + x² + 2u - 0.2xu + x²u
%       y = 2x
%
%   which has an analytic solution to the past and future energy functions.
%
%   Excluding the higher order g term gives the same system as in [1]. The polynomial
%   version is used in [2].
%
%   Reference: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki,
%               “Scalable computation of energy functions for nonlinear
%               balanced truncation,” Computer Methods in Applied Mechanics
%               and Engineering, vol. 427, p. 117011, Jul. 2024, doi:
%               10.1016/j.cma.2024.117011
%              [2] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗∞
%               energy functions for polynomial control-affine systems,"
%               IEEE Transactions on Automatic Control, pp. 1–13, 2024,
%               doi: 10.1109/tac.2024.3494472
%
%   Part of the NLbalancing repository.
%%

A = -2;
B = 2;
C = 2;
N = 1;
G1 = -0.2;
G2 = 0.2;

f = {A, N};
g = {B, G1, G2};
h = {C};

end
