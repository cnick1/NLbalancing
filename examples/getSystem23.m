function [f, g, h] = getSystem23(linear)
%getSystem23  Returns the 3D model from Clancy Rowley's talk. Need to find original sources.
%
%   Usage:  [f,g,h] = getSystem23(linear)
%
%   Inputs:
%       linear - boolean, whether to use the linear model or nonlinear
%                model from the talk
%
%   References: [1]
%
%%

if nargin < 1 
    linear = false;
end


n = 3; 
x = sym('x', [1, n]).'; syms(x);

if linear
    fsym = [-x1 + 100*x3;
        -2*x2 + 100*x3;
        -5*x3];
    gsym = [1;1;1];
    hsym = x1+x2+x3;
else
    fsym = [-x1 + 20*x1*x3;
        -2*x2 + 20*x2*x3;
        -5*x3];
    gsym = [1;1;1];
    hsym = x1+x2+x3;
end


[f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, 3);

end
