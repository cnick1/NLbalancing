function [f, g, h] = getSystem10()
%getSystem10  Generates a 1D polynomial system with transcritical bifurcation.
%
%   Usage:  [f,g,h] = getSystem10()
%
%   Outputs:    f,g,h  - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The dynamics correspond to the input-output system
%       ẋ = 1/10 x - 2 x² + x³ + u
%       y = x
%
%   The state equation is similar to the model from [1].
%
%   References: [1] J. Borggaard and L. Zietsman, “Computation of nonlinear feedback
%               for flow control problems,” in 2018 American Control Conference
%               (ACC), Jun. 2018, doi: 10.23919/acc.2018.8431410.
%
%   Part of the PPR repository.
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
