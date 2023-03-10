function [A, B, C, N, Q, f, g, h] = getSystem6()
%getSystem6  Generates a unicycle quadratic-bilinear system for testing energy functions.
%
%   Usage:  [A,B,C,N,Q] = getSystem6()
%        or [~,~,~,~,~,f,g,h] = getSystem6()
%
%   The "matrices" correspond to the quadratic bilinear input-output system
%
%       \dot(x1) = u1 cos x3
%       \dot(x2) = u1 sin x3  - x2
%       \dot(x3) = u2
%             y = x
%
%   upon 'linearizing' f(x) and g(x) (leading to a bilinear system rather
%   than linear, in which case g(0) should be taken). The
%   cell arrays f, g, and h can also be used as outputs. 
%
%   Part of the NLbalancing repository.
%%

A = [0 0 0; 0 -1 0; 0 0 0];
B = [1 0; 0 0; 0 1];
C = speye(3);
N = sparse(3, 3);
Q = [0 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 0];

f = {A,N}; 
g = {B,Q}; 
h = {C}

end
