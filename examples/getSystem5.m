function [f, g, h] = getSystem5()
%getSystem5  Generates a unicycle quadratic-bilinear system for testing energy functions.
%
%   Usage:  [f,g,h] = getSystem5()
%
%   The "matrices" correspond to the input-output system
%
%       \dot{x1} = u1 cos x3
%       \dot{x2} = u1 sin x3  - x2
%       \dot{x3} = u2
%             y = x
%
%   upon polynomial approximation to f(x) and g(x). The cell arrays f, g, and h
%   can also be used as outputs.
%
%   Part of the NLbalancing repository.
%%

A = [0 0 0;
     0 -1 0;
     0 0 0];
B = [1 0;
     0 0;
     0 1];
C = speye(3);
N = sparse(3, 3);
G = [0 0 0 0 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 0];

f = {A, N};
g = {B, G};
h = {C};

end
