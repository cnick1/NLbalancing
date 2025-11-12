function [f, g, h, FofXU] = getSystem31(degree)
%getSystem31 Generates a polynomial 2D stable pendulum model
%
%   Usage:  [f,g,h,FofXU] = getSystem31()
%
%   Inputs:    degree - desired degree of the polynomial approximation
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%                FofXU - Function handle for the true dynamics ẋ = f(x,u)
%
%   Description: We consider the pendulum as described in [1-3] by
%
%             ẍ + G/L sin(x) + k / m L² x + b/mL² ẋ = 1/mL² u(t)
%
%   This can be put in first-order form with x₁ = x, x₂ = ẋ as
%
%             ẋ₁ = x₂
%      (1)    ẋ₂ = -G/L sin(x₁) - k / m L² x₁ - b/mL² x₂ + u(t)
%              y = x₁
%
%   See Ex. 26 for a possible alternative state space system.
%
%   References: [1] A. J. Newman and P. S. Krishnaprasad, “Computation for
%                   nonlinear balancing,” University of Maryland, College
%                   Park, 1998.
%               [2] A. J. Newman and P. S. Krishnaprasad, “Computing
%                   balanced realizations for nonlinear systems,” University
%                   of Maryland, College Park, 2000.
%               [3] A. J. Newman, “Modeling and reduction with applications
%                   to semiconductor processing,” University of Maryland,
%                   College Park, 1999.
%
%%
if nargin < 1
    degree = 9;
end

n = 2;

%       ẋ₁ = x₂
%  (1)  ẋ₂ = -G/L sin(x₁) - k / m L^2 x₁ - b/mL^2 x₂ + u(t)
%        y = x₁
g=10; L=20; m=1/40; b=2; k=1;

FofXU = @(x,u) [x(2); -g/L*sin(x(1)) - k/(m*L^2)*x(1)  - b/(m*L^2)*x(2) + u/(m*L^2)];

A = [0 1;
    -g/L-k/(m*L^2) -b/(m*L^2)];
F2 = sparse(n,n^2);
f = {A, F2};
for i = 1:(degree - 1) / 2
    f{end + 1} = sparse(2, 2 ^ (2 * i + 1));
    f{end}(2, 1) = (-1) ^ i / factorial(2 * i + 1) * (-g/L);
    f{end + 1} = sparse(2, 2 ^ (2 * i + 2));
end
f = f(1:degree);

B = [0; 1/(m*L^2)];
C = [1 0];

g = {B};
h = {C};

end
