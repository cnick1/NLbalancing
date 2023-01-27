function [A,B,C,N,Q] = getSystem5()
%getSystem5  Generates a simple quadratic system for testing energy functions.
%
%   Usage:  [A,B,C,N,Q] = getSystem5()
%
%   The "matrices" correspond to the quadratic bilinear input-output system
%
%       \dot(x) = -2x + x^2 + (2+2x)u
%             y = 2x
%
%   for which there is an analytic solution to the past and future energy
%   functions.
%
%   Reference: Nonlinear Balanced Truncation Model Reduction: 
%        Part 1-Computing Energy Functions, by Kromer, Gugercin, and Borggaard.
%        arXiv.
%
%   Part of the NLbalancing repository.
%%

    A = -2;
    B =  2;
    C =  2;
    N =  1;
    Q =  2;

end
