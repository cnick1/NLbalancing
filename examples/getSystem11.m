function [f, g, h] = getSystem11(degree, m, L)
%getSystem11  Generates a 2D inverted pendulum model with "degree"
%             polynomial sin approximation.
%
%   Usage:  [f,g,h] = getSystem11(degree, m, L)
%
%   Inputs:    degree - desired degree of the polynomial approximation
%                   m - mass
%                   L - pendulum length
%
%   Outputs:    f,g,h  - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%               FofXU  - function handle for the function
%                               ẋ = f(x,u)
%                                 = f(x) + g(x) u
%                        Typically you could construct this with
%                        kronPolyEval(f,x), etc., but since the model here
%                        is non-polynomial, it is more correct to use the
%                        non-polynomial model, e.g. for simulations.
%
% Description: This example is for a pendulum made from a uniform rod. For
% the idealized pendulum that is non-dimensionalized, use getSystem26
% instead. The polynomial approximation to the pendulum is
%           ẋ₁ = x₂
%           ẋ₂ = 3u/(mL²) + 3g/(2L) ( x₁ - x₁³/6 + x₁⁵/120 - x₁⁷/5040 + x₁⁹/362880 + ... )
%           y  = x₁
%
%   Reference: [1] N. A. Corbin and B. Kramer, “Computing solutions to the
%               polynomial-polynomial regulator problem,” in 2024 63rd IEEE
%               Conference on Decision and Control, Dec. 2024
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


F = @(x) [x(2); 3 * gravity / (2 * L) * sin(x(1))];
G = @(x) B;
FofXU = @(x,u) (F(x) + G(x)*u);


end
