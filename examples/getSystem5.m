function [f, g, h] = getSystem5(eps, degree)
%getSystem5  Generates a unicycle quadratic-bilinear system for testing energy functions.
%
%   Usage:  [f,g,h] = getSystem5()

%   Inputs:
%       eps     -  coefficient on the stabilizing 'wind' in the left-right direction
%       degree  -  degree d taylor approximation to the polynomial input g(x)
%
%   Outputs:
%       f,g,h   - cell arrays containing the polynomial coefficients
%                 for the drift, input, and output.
%
%   Background: The outputs correspond to f(x), g(x), and h(x) for the input-output system
%
%       \dot(x1) = u1 cos x3
%       \dot(x2) = u1 sin x3  - eps*x2
%       \dot(x3) = u2
%             y = x
%
%   Part of the NLbalancing repository.
%%

if nargin < 2
    degree = 7;
    if nargin < 1
        eps = 0.1;
    end
end

A = [0 0 0;
     0 -eps 0;
     0 0 0];

B = [1 0;
     0 0;
     0 1];

C = eye(3);

F2 = sparse(3, 9);

% Form permuted Gs since it is easier, then permute
g = {B};
for i = 0:(degree - 1) / 2
    % sin(x3)
    g{end + 1} = sparse(3, 3 ^ (2 * i + 1)); % Only form first half so that I can index the end, then append zeros afterward
    g{end}(2, end) = (-1) ^ i / factorial(2 * i + 1);
    g{end} = [g{end}, sparse(3, 3 ^ (2 * i + 1))];
    % Permute
    g{end} = g{end} * perfectShuffle(3 ^ (2 * i + 1), 2);

    % cos(x3)
    g{end + 1} = sparse(3, 3 ^ (2 * i + 2));
    g{end}(1, end) = (-1) ^ i / factorial(2 * i + 2);
    g{end} = [g{end}, sparse(3, 3 ^ (2 * i + 2))];
    % Permute
    g{end} = g{end} * perfectShuffle(3 ^ (2 * i + 2), 2);
end
g = g(1:degree);

f = {A, F2};
h = {C};

end
