function [A, B, C, N, Q, f, g, h] = getSystem5()
%getSystem5  Generates a simple quadratic-bilinear system for testing energy functions.
%
%   Usage:  [A,B,C,N,Q] = getSystem5()
%        or [~,~,~,~,~,f,g,h] = getSystem5()
%
%   The "matrices" correspond to the quadratic-bilinear input-output system
%
%       \dot(x) = -2x + x^2 + (2+2x)u
%             y = 2x
%
%   for which there is an analytic solution to the past and future energy functions.
%
%   Excluding the bilinear term gives the same system as getSystem1(). The
%   cell arrays f, g, and h can also be used as outputs. 
%
%   Reference: Nonlinear Balanced Truncation Model Reduction:
%        Part 1-Computing Energy Functions, by Kromer, Gugercin, and Borggaard.
%        arXiv.
%
%   Part of the NLbalancing repository.
%%

A = -2;
B = 2;
C = 2;
N = 1;
Q = 2;

f = {A,N}; 
g = {B,Q}; 
h = {C}

end
