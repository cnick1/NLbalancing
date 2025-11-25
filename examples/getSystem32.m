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
%   Consider the following coordinate transformation, which we apply to the
%   balanced realization of the model:
%
%           x = 𝟁(z) = [z₁;   z₂;   z₃ + z₁² + z₂² + z₁³]
%
%   The dynamics in the z coordinates are then
%
%   ż₁ = -0.739 z₁ + 1.57 z₂ - 0.172 z₃ - 0.172 z₁² - 0.172 z₂² - 0.172 z₁³ + 5.09 u
%   ż₂ = -1.57 z₁ - 6.26 z₂ + 1.72 z₃ + 1.72 z₁² + 1.72 z₂² + 1.72 z₁³ + 4.82 u
%   ż₃ = -0.172 z₁ - 1.72 z₂ - 1.0 z₃ + 0.343 z₁z₃ - 3.43 z₂z₃ + 0.476 z₁² + 11.5 z₂² + 0.343 z₁z₂² - 8.13 z₁²z₂ + 0.515 z₁²z₃ + 1.56 z₁³ - 3.43 z₂³ + 0.515 z₁²z₂² - 3.43 z₁³z₂ + 0.859 z₁⁴ + 0.515 z₁⁵ + (0.598 - 15.3 z₁² - 10.2 z₁ - 9.64 z₂) u
%    y = 5.09 z₁ - 4.82 z₂ +0.598 z₃ + 0.597 z₁³ + 0.597 z₁² + 0.597 z₂²
%
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
    nvp.transform = true
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
    % will do here. To start, let's compute that balanced realization.
    [fbal,gbal,hbal,Tbal] = getBalanceThenReduceRealization(f,g,h,eta=0,degree=1);
    
    % Now we will transform by the coordinate transformation
    %        x = 𝟁(z) = [z₁;   z₂;   z₃ + z₁² + z₂² + z₁³]
    
    a = 1; z = sym('z', [1, n]).';
    [Tnl,~,~] = approxPolynomialDynamics([z(1); z(2); z(3) + a*(z(1)^2 + z(2)^2) + a*(z(1)^3)], gsym, z(1), z, 3);
    
    [f,g,h] = transformDynamics(fbal,gbal,hbal,Tnl,degree=5);
    % [f,g,h] = transformDynamics(f,g,h,transformationInverse(Tbal),degree=5); % Can also reverse the linear transformation
    
end

end
