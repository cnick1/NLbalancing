function [A, B, C, N, f, g, h] = getSystem1()
%getSystem1  Generates a simple 1D polynomial system for testing energy functions.
%
%   Usage:  [A,B,C,N] = getSystem1()
%        or [A,B,C,N,f,g,h] = getSystem1()
%
%   The "matrices" correspond to the input-output system
%
%       \dot(x) = -2x + x^2 + 2u - 0.2xu + x^2u
%             y = 2x
%
%   for which there is an analytic solution to the past and future energy functions.
%
%   Excluding the higher order g term gives the same system as in [1]. The
%   cell arrays f, g, and h can also be used as outputs.
%
%   Reference: [1] Nonlinear Balanced Truncation Model Reduction:
%        Part 1-Computing Energy Functions, by Kromer, Gugercin, and Borggaard.
%        arXiv.
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
