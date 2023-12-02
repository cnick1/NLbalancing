function [f, g, h] = getSystem15(degree)
%getSystem15  Returns a polynomial approximation to the 4D model of a double pendulum.
%
%   Usage:  [f,g,h] = getSystem15(degree)
%
%   Inputs:
%       degree - desired degree of the computed polynomial dynamics.
%
%   Outputs:
%       f,g,h  - Cell arrays containing the polynomial coefficients for the
%                drift, input, and output (generalizations containing A,B,C)
%
%   Background: This model has been used several times in the literature.
%       Let the 2 x 2 mass matrix be given by the entries
%           m11       = m1 l1^2 + m2 l1^2 + m2 l2^2 + 2 m2 l1 l2 cos x2
%           m12 = m21 = m2 l2^2 + m2 l1 l2 cos x2
%           m22       = m2 l2^2
%       The mass matrix and its inverse are then
%           M(x) = [m11, m12;    M^(-1)(x) = _______1̲_______  [m22, -m21;
%                   m21, m22]               (m11m22 - m12m21) -m12,  m11]
%
%       The potential energy of the system is
%           V(x) = - m1 g l1 cos x1 - m2 g (l1 cos x1 + l2 cos(x1 + x2))
%
%       and the kinetic energy is (in [1,2], they have a typo  missing the 1/2)
%           T(ẋ) = 1/2 ẋ.' * M * ẋ  
%
%       The Lagrangian is L(x,ẋ) = T(ẋ) - V(x). The Euler-Lagrange equations
%       are
%           d/dt(∂L/∂ẋ) - ∂L/∂x = Q

%       Here, our generalized forces are the control input and damping, given
%       by Q = [u - μ1ẋ1; - μ2ẋ2]. For our problem, ∂L/∂ẋ = ∂T/∂ẋ = Mẋ, so
%       d/dt(∂L/∂ẋ) = d/dt(Mẋ); since the mass matrix is state dependent, this
%       is d/dt(Mẋ) = Mẍ + Mdot ẋ. So the Euler-Lagrange equations become
%           Mẍ + Ṁ ẋ - ∂L/∂x = Q
%
%       Converting the second order system to first order, we define our state
%       as x = [x1 x2 ẋ1 ẋ2], so ẋ = [ẋ1 ẋ2 ẍ1 ẍ2]. Thus the governing 
%       equations are
%           ẋ = [                  x3; 
%                                   x4;
%                 M^(-1)( ∂L/∂(x1 x2) - Mdot [x3; x4] + Q)]
%
%       The output is taken as the horizontal and vertical positions of the
%       second mass. These are
%           y = [l1 * sin(x1) + l2 * sin(x1 + x2);
%                l1 * cos(x1) + l2 * cos(x1 + x2)];
%
%       NOTE: This is the model as presented in [1,2]. However, the damping here
%       is not physically correct, because it neglects to apply the appropriate
%       damping resistance to both bars (equal and opposite). Note how μ1 only
%       enters the first equation, and μ2 only enters the second equation.
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
