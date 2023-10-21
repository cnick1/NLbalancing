function [f, g, h] = getSystem11(degree, m, L)
%getSystem11  Generates a 2D inverted pendulum model with "degree"
%             polynomial sin approximation.
%
%   Usage:  [f,g,h] = getSystem11(degree, m, L)
%
%              \dot{x}_1 & = x_2                                                                                                                     \\
%              \dot{x}_2 & = \frac{3T}{m\L^2} + \frac{3g}{2\L} (x_1-\frac{x_1^3}{6}+\frac{x_1^5}{120}-\frac{x_1^7}{5040} + \frac{x_1^9}{362880})
%                   y    & = x_1
%
%   Reference: [1] N. A. Corbin and B. Kramer, “The polynomial-polynomial regulator:
%              computing feedback for polynomially nonlinear systems with polynomial
%              performance indexes,” 2023.
%
%%

if nargin < 3
    if nargin < 2
        if nargin < 1
            degree = 7;
        end
        m = 1;
    end
    L = 1;
end

gravity = 9.81;

A = [0, 1; 3 * gravity / (2 * L), 0];
F2 = sparse(2, 2 ^ 2);

f = {A, F2};
for i = 1:(degree - 1) / 2
    f{end + 1} = sparse(2, 2 ^ (2 * i + 1));
    f{end}(2, 1) = 3 * gravity / (2 * L) * (-1) ^ i / factorial(2 * i + 1);
    f{end + 1} = sparse(2, 2 ^ (2 * i + 2));
end
f = f(1:degree);

B = [0; 3 / (m * L ^ 2)]; % Torque control at pivot point

C = [1, 0]; % Measure angle of pendulum

g = {B};
h = {C};

end
