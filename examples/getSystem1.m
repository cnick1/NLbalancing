function [f, g, h] = getSystem1()
%getSystem1  Generates a simple 1D polynomial system for testing energy functions.
%
%   Usage:  [f, g, h] = getSystem1()
%
%   The dynamics correspond to the input-output system
%
%       \dot{x} = -2x + x^2 + 2u - 0.2xu + x^2u
%             y = 2x
%
%   which has an analytic solution to the past and future energy functions.
%
%   Excluding the higher order g term gives the same system as in [1]. The polynomial
%   version is used in [2].
%
%   Reference: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki, “Nonlinear
%               balanced truncation: Part 1—computing energy functions,” arXiv,
%               Dec. 2022. doi: 10.48550/ARXIV.2209.07645
%              [2] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” 2023.
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
