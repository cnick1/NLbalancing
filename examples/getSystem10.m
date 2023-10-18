function [A, B, C, F2, f, g, h] = getSystem10()
%getSystem10  Generates a simple 1D polynomial system for testing basin of attraction.
%
%   Usage:  [A,B,C,F2] = getSystem10()
%        or [A,B,C,F2,f,g,h] = getSystem10()
%        or [~,~,~,~,f,g,h] = getSystem10()
%
%   The "matrices" correspond to the input-output system
%
%       \dot(x) = -1/10 x - 4 x^2 + x^3 + u
%             y = x
%
%   for which there is an analytic solution to the past and future energy functions.
%
%   The state equation is from [1].
%
%   References: [1] J. Borggaard and L. Zietsman, “Computation of nonlinear feedback
%   for flow control problems,” in 2018 American Control Conference (ACC), Jun. 2018,
%   doi: 10.23919/acc.2018.8431410.
%
%   Part of the NLbalancing repository.
%%

A = 1/10;
F2 = -2;
F3 = 1;
B = 1;
C = 1;

f = {A, F2, F3};
g = {B};
h = {C};

end
