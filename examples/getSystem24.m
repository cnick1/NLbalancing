function [f, g, h] = getSystem24(linear)
%getSystem24 Returns the 2D nonlinear model based on [1,2]. 
%
%   Usage:  [f,g,h] = getSystem24()
%
%   The dynamics correspond to the input-output system
%
%       xdot_1 = −2 x1 + 20 x1 x2 + u,
%       xdot_2 = −5 x2 + u,
%            y = x1 + x2
%
%   If the option linear is enabled, the model from [3, Section 5.6.1] is 
%   returned instead. This model replaces the nonlinear interaction between 
%   x3 and x1 & x2 with a linear interaction. The model from [3] inspired 
%   the model in [1,2]. 
% 
%       xdot_1 = −2 x1 + 100 x2 + u,
%       xdot_2 = −5 x2 + u,
%            y = x1 + x2,
%
%   Inputs:     linear - boolean, whether to use the linear model from [3]  
%                        or nonlinear model from [1,2]
% 
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output 
%
%   References: 
%
%%

if nargin < 1 
    linear = false;
end


n = 2; 
x = sym('x', [1, n]).'; syms(x);

if linear
    fsym = [-2*x1 + 100*x2;
        -5*x2];
    gsym = [1;1];
    hsym = x1+x2;
else
    fsym = [-2*x1 + 20*x1*x2;
        -5*x2];
    gsym = [1;1];
    hsym = x1+x2;
end


[f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, 3);

end
