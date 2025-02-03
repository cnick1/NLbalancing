function [f, g, h, fsym, gsym, hsym] = getSystem15(degree)
%getSystem15  Returns a polynomial approximation to the 4D model of a double pendulum.
%
%   Usage:  [f,g,h] = getSystem15(degree)
%
%   Inputs:
%       degree - desired degree of the computed polynomial dynamics.
%
%   Outputs:    f,g,h  - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output (generalizations
%                        containing A,B,C)
%       fsym,gsym,hsym - symbolic handles for the true dynamics
%
%   Description: This model has been used several times in the literature.
%   The state variables are
%       x₁ - angular position of mass 1, relative to vertical
%       x₂ - angular position of mass 2, relative to the angle of mass 1
%            (straight in line would be x₂=0)
%       ẋ₁ - angular velocity of mass 1, relative to vertical
%       ẋ₂ - angular velocity of mass 2, relative to the angle of mass 1
%   Let the 2 x 2 mass matrix be given by the entries
%       m₁₁       = m₁ l₁² + m₂ l₁² + m₂ l₂² + 2 m₂ l₁ l₂ cos x₂
%       m₁₂ = m₂₁ = m₂ l₂² + m₂ l₁ l₂ cos x₂
%       m₂₂       = m₂ l₂²
%   The mass matrix and its inverse are then
%       M(x) = [m₁₁, m₁₂;    M⁻¹(x) = _______1_______  [m₂₂, -m₂₁;
%               m₂₁, m₂₂]            (m₁₁m₂₂ - m₁₂m₂₁) -m₁₂,  m₁₁]
%
%   The potential energy of the system is
%       V(x) = - m₁ g l₁ cos x₁ - m₂ g (l₁ cos x₁ + l₂ cos(x₁ + x₂))
%
%   and the kinetic energy is (in [1,2], they have a typo missing the 1/2)
%       T(ẋ) = 1/2 ẋ.' * M * ẋ
%
%   The Lagrangian is L(x,ẋ) = T(ẋ) - V(x). The Euler-Lagrange equations are
%       d/dt(∂L/∂ẋ) - ∂L/∂x = Q
%
%   Here, our generalized forces are the control input and damping, given
%   by Q = [u - μ₁ẋ₁; - μ₂ẋ₂]. For our problem, ∂L/∂ẋ = ∂T/∂ẋ = Mẋ, so
%   d/dt(∂L/∂ẋ) = d/dt(Mẋ); since the mass matrix is state dependent, this
%   is d/dt(Mẋ) = Mẍ + Ṁẋ. So the Euler-Lagrange equations become
%       Mẍ + Ṁ ẋ - ∂L/∂x = Q
%
%   Converting the second order system to first order, we define our state
%   as x = [x₁ x₂ ẋ₁ ẋ₂], so ẋ = [ẋ₁ ẋ₂ ẍ₁ ẍ2]. Thus the governing
%   equations are
%       ẋ = [                  x₃;
%                              x₄;
%              M⁻¹( ∂L/∂(x₁x₂) - Ṁ [x₃; x₄] + Q) ]
%
%   The output is taken as the horizontal and vertical positions of the
%   second mass. These are
%       y = [l₁ * sin(x₁) + l₂ * sin(x₁ + x₂);
%            l₁ * cos(x₁) + l₂ * cos(x₁ + x₂)]
%
%   References: [1] K. Fujimoto and D. Tsubakino, “On computation of
%               nonlinear balanced realization and model reduction,” in
%               2006 American Control Conference, IEEE, 2006. doi:
%               10.1109/acc.2006.1655399
%               [2] K. Fujimoto and D. Tsubakino, “Computation of nonlinear
%               balanced realization and model reduction based on Taylor
%               series expansion,” Systems & Control Letters, vol. 57, no.
%               4, pp. 283–289, Apr. 2008, doi: 10.1016/j.sysconle.2007.08.015
%
%
%
%%

if nargin < 1
    degree = 3;
end

n = 4;
x = sym('x', [1, n]).';
syms(x);

G = 9.8;
m1 = 1; m2 = 1;
l1 = 1; l2 = 1;
mu1 = 1; mu2 = 1;

% define joint 2 angle relative to joint 1 angle, as in [1,2]
V =- m1 * G * l1 * cos(x1) ... % m1 potential energy
    - m2 * G * (l1 * cos(x1) + l2 * cos(x1 + x2)); % m2 potential energy

% Using velocities as x3 and x4
m11 = m1 * l1 ^ 2 + m2 * l1 ^ 2 + m2 * l2 ^ 2 + 2 * m2 * l1 * l2 * cos(x2);
m12 = m2 * l2 ^ 2 + m2 * l1 * l2 * cos(x2);
m22 = m2 * l2 ^ 2;

M = [m11, m12;
    m12, m22];
Minv = 1 / (m11 * m22 - m12 ^ 2) * ...
    [m22, -m12;
    -m12, m11];
Mdot =- [2 * m2 * l1 * l2 * sin(x2) * x4, m2 * l1 * l2 * sin(x2) * x4 ;
    m2 * l1 * l2 * sin(x2) * x4 , 0];

T = 0.5*[x3; x4].' * M * [x3; x4];
% T = [x3; x4].' * M * [x3; x4];

fsym = [x3;
    x4;
    Minv * (gradient(T - V, [x1; x2]) - Mdot * [x3; x4] - [mu1 * x3; mu2 * x4])];

gsym = [0; 0; Minv * [1; 0]];
hsym = [l1 * sin(x1) + l2 * sin(x1 + x2);
    l1 * (1-cos(x1)) + l2 * (1-cos(x1 + x2))];

[f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, degree);

if false
    %% port Hamiltonian approach
    
    q = [x1; x2]; qdot = [x3; x4];
    
    L = T - V;
    
    % Generalized momenta conjugate to our choice of generalized coordinates
    % simplify(gradient(L, qdot) - M*qdot) % Check that M satisfies p = M*qdot
    
    % Define Hamiltonian using Legendre transform
    % simplify(T+V - (p.' * qdot - L)) % Check that H=T+V
    
    %% Port-Hamiltonian approach to deriving equations of motion
    % redefine Hamiltonian state
    clear x3 x4
    syms p1 p2
    y = [q; p1; p2];
    
    T = 1/2 * [p1; p2].' * Minv * [p1; p2];
    
    H = T+V;
    D = [mu1, 0;
        0, mu2];
    
    J = [zeros(n / 2), eye(n / 2);
        -eye(n / 2), zeros(n / 2)];
    R = [zeros(n / 2), zeros(n / 2);
        zeros(n / 2), D];
    
    fsym = (J - R) * gradient(H, y);
    gsym = [0; 0; 1; 0];
    hsym = [l1 * sin(x1) + l2 * sin(x1 + x2);
        l1 * cos(x1) + l2 * cos(x1 + x2)];
    
    [f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, y, degree);
end
end
