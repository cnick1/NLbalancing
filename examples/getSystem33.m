function [f, g, h] = getSystem33()
%getSystem32 Returns Hirsch's repressilogilator
%
%   Usage:  [f,g,h] = getSystem33()
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
%   References: [1] 
%
%%
n = 5; x = sym('y', [1, n]).';
K1 = 0.25; K2 = 0.25; K3 = 0.25; K4 = 0.25; K5 = 0.25;
n1 = 4.0; n2 = 4.0; n3 = 4.0; n4 = 2.0; n5 = 2.0;
fsymtemp = [-x(1) + 1/(1+(x(3)/K3)^n3);
    -x(2) + 1/(1+(x(1)/K1)^n1);
    -x(3) + 1/(1+(x(2)/K2)^n2);
    -x(4) + 1/(1+(x(5)/K5)^n5);
    -x(5) + 1/(1+(x(4)/K4)^n4);
];
gsymtemp = [0;0;0;1/(1+(x(2)/K2)^n2);1/(1+(x(3)/K3)^n3)];
hsymtemp = x;

syms x1 x2 x3 x4 x5 a b 
y1 = x1 + a;
y2 = x2 + a;
y3 = x3 + a;
y4 = x4 + b;
y5 = x5 + b;
y = [y1;y2;y3;y4;y5];

fsym = subs(fsymtemp,x,y);
gsym = subs(gsymtemp,x,y);
hsym = subs(hsymtemp,x,y);
x = sym('x', [1, n]).';

eqs = subs(fsym,x,[0;0;0;0;0]);
sol1 = double(solve(eqs(1) == 0, a)); 
sol2 = double(solve(eqs(5) == 0, b)); 

fsym = subs(fsym,{a, b},{sol1(1), sol2(1)});
gsym = subs(gsym,{a, b},{sol1(1), sol2(1)});
hsym = x;

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
