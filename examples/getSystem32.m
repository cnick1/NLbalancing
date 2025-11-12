function [f, g, h] = getSystem32(nvp)
%getSystem32  Returns Holmes' 3D linear balanced truncation model.
%
%   Usage:  [f,g,h] = getSystem32()
%           [f,g,h] = getSystem32(transform=true)
%
%   Name/Value Pair Options:     
%   transform - If set to true, a nonlinear transformation is applied to
%               the linear model to make it nonlinear. The nonlinear
%               balancing algorithm should be able to "undo" the
%               transformation, as the balanced realization is linear. 
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The dynamics correspond to the input-output system
%
%           ẋ₁ = −x₁ + 100 x₃ + u,
%           ẋ₂ = −2 x₂ + 100 x₃ + u,
%           ẋ₃ = −5 x₃ + u,
%            y = x₁ + x₂ + x₃.
%
%   This model was derived in [1, Section 5.6.1]] as an example of a reducable
%   model. The output contains all of the states, and the input affects all of 
%   the states. The third state is not exactly decoupled from the others, and
%   while it decays quickly, it also drives the other states very strongly.
%   Nonetheless, there is a two-dimensional model that captures the
%   input-output behavior (transfer function, impulse and step response,
%   etc.) very well, and it is not just the first two states. 
%
%   As a linear model, this example can be used to validate that the
%   nonlinear balanicng algorithms reduce to the established linear methods
%   when applied to linear systems. However, we can also apply a nonlinear
%   transformation to this model in order to produce a nonlinear model. 
%   The balanced realization should still be the same linear balanced
%   realization, so the nonlinear balancing algorithm should in effect
%   "undo" the nonlinear transformation we apply. 
%
%   References: [1] P. Holmes, J. L. Lumley, G. Berkooz, and C. W. Rowley,
%                   Turbulence, coherent structures, dynamical systems and
%                   symmetry. Cambridge University Press, 2012. doi:
%                   10.1017/cbo9780511919701.
%
%   See also: getSystem23
%%
arguments
    nvp.transform = false
end

n = 3; x = sym('x', [1, n]).'; 

fsym = [-x(1) + 100*x(3);
    -2*x(2) + 100*x(3);
    -5*x(3)];
gsym = [1;1;1];
hsym = x(1)+x(2)+x(3);

[f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, 3);

if nvp.transform
% The balanced realization is a 2D plane in a 3D space for this linear
% system. Imagine curving that 2D plane into a manifold; that is what we
% will do here. To start, let's compute that 2D balanced subspace 
L = lyapchol(f{1},g{1}); R = lyapchol(f{1}.',h{1}.');
[U,S,V] = svd(L*R.');
T = {invertibleMatrix(L.' * U / S.^(1/2), S.^(1/2) \ V.' * R)};


[ft, gt, ht] = transformDynamics(f, g, h, T);

end

end
