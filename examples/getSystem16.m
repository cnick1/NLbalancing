function [f, g, h] = getSystem16(degree)
%getSystem16  Returns a polynomial approximation to the 6D model of a triple pendulum.
%
%   Usage:  [f,g,h] = getSystem16(degree)
%
%   Inputs:
%       degree - desired degree of the computed polynomial dynamics.
%
%   Outputs:
%       f,g,h  - Cell arrays containing the polynomial coefficients for the
%                drift, input, and output (generalizations containing A,B,C)
%
%   Background: We consider a triple pendulum with rigid uniform slender
%       links. Each joint also has an associated spring stiffness and
%       damping. We derive the equations of motion using either a
%       Lagrangian or Hamiltonian approach. We begin by defining the
%       kinetic end potential energies of the system. The kinetic energy is
%       composed of
%           - m1 kinetic energy (fixed point)
%           - m2 rotational kinetic energy
%           - m2 translational kinetic energy
%           - m3 rotational kinetic energy
%           - m3 translational kinetic energy
%           (see the code for the specific formulas; essentially one needs
%           to use vector differentiation/geometry, along with Ig=m l^2/12)
%       The potential energy is composed of
%           - m1 potential energy (gravity)
%           - k1 potential energy (spring)
%           - m2 potential energy (gravity)
%           - k2 potential energy (spring)
%           - m3 potential energy (gravity)
%           - k3 potential energy (spring)
%           (see the code for the specific formulas; essentially one needs
%           to use geometry to find the position of each center of gravity,
%           and for the springs E = 1/2 k Δx^2)
%
%       With the kinetic energy T and potential energy V, we can either use
%       Lagrange's equations or a port-Hamiltonian approach
%           d/dt(∂L/∂ẋ) - ∂L/∂x = Q     or      ẋ = (J-R)∂H/∂x + Bu
%       where the former requires inverting and differentiating a mass
%       matrix to convert to first-order (state-space) equations.
%
%   References: [1] A. J. Krener, “Reduced order modeling of nonlinear
%               control systems,” in Analysis and Design of Nonlinear
%               Control Systems, Springer Berlin Heidelberg, 2008, pp.
%               41–62. doi: 10.1007/978-3-540-74358-3_4.
%
%
%%

if nargin < 1
    degree = 3;
end

n = 6;
x = sym('x', [1, n]).';
syms(x);

G = .5;

m1 = 1; m2 = 1; m3 = 1;
l1 = 2; l2 = 2; l3 = 2;
k1 = 3; k2 = 3; k3 = 3;
mu1 = 0.5; mu2 = 0.5; mu3 = 0.5;

% define positions from fixed vertical, see e.g. https://scienceworld.wolfram.com/physics/DoublePendulum.html
% Using velocities as x4 x5 x6
T = 0.5 * 1/3 * m1 * l1 ^ 2 * x4 ^ 2 ... % m1 kinetic energy (fixed point)
    + 0.5 * 1/12 * m2 * l2 ^ 2 * x5 ^ 2 ... % m2 rotational kinetic energy
    + 0.5 * m2 * (l1 ^ 2 * x4 ^ 2 +1/4 * l2 ^ 2 * x5 ^ 2 + l1 * l2 * x4 * x5 * cos(x2 - x1)) ... % m2 translational kinetic energy
    + 0.5 * 1/12 * m3 * l3 ^ 2 * x6 ^ 2 ... % m3 rotational kinetic energy
    + 0.5 * m3 * (l1 ^ 2 * x4 ^ 2 + l2 ^ 2 * x5 ^ 2 +1/4 * l3 ^ 2 * x6 ^ 2 ... % m3 translational kinetic energy
    + 2 * l1 * l2 * x4 * x5 * cos(x2 - x1) + l1 * l3 * x4 * x6 * cos(x3 - x1) ...
    + l2 * l3 * x5 * x6 * cos(x3 - x2));

V =- m1 * G * l1 / 2 * cos(x1) ... % m1 potential energy (gravity)
    +1/2 * k1 * x1 ^ 2 ... % k1 potential energy (spring)
    - m2 * G * (l1 * cos(x1) + l2 / 2 * cos(x2)) ... % m2 potential energy (gravity)
    +1/2 * k2 * (x2 - x1) ^ 2 ... % k2 potential energy (spring)
    - m3 * G * (l1 * cos(x1) + l2 * cos(x2) + l3 / 2 * cos(x3)) ... % m3 potential energy (gravity)
    +1/2 * k3 * (x3 - x2) ^ 2; % k3 potential energy (spring)

q = [x1; x2; x3]; qdot = [x4; x5; x6];

L = T - V;

% Generalized momenta conjugate to our choice of generalized coordinates
p = gradient(L, qdot);

% Define Hamiltonian using Legendre transform
H = p.' * qdot - L;

%         simplify(T + V - H) % Check: Hamiltonian turns out to be total energy

conversionMethod = 1;
switch conversionMethod
    case 1
        %% 2nd order mechanical system method
        % Derive a mass matrix such that the momentum is p = M*qdot
        % Hard coded!
        M = [28/3, 6 * cos(x1 - x2), 2 * cos(x1 - x3);
             6 * cos(x1 - x2), 16/3, 2 * cos(x2 - x3);
             2 * cos(x1 - x3), 2 * cos(x2 - x3), 4/3];
        %         simplify(p - M*qdot) % Check that M satisfies p = M*qdot
        %                 Minv = inv(M); % yikes lol
        Mdot =- [0, 6 * sin(x1 - x2) * (x4 - x5), 2 * sin(x1 - x3) * (x4 - x6);
                 6 * sin(x1 - x2) * (x4 - x5), 0, 2 * sin(x2 - x3) * (x5 - x6);
                 2 * sin(x1 - x3) * (x4 - x6), 2 * sin(x2 - x3) * (x5 - x6), 0];

        fsym = [qdot;
                M \ (gradient(L, q) - Mdot * qdot - [(mu1 + mu2) * x4 - mu2 * x5; -mu2 * x4 + (mu2 + mu3) * x5 - mu3 * x6; - mu3 * x5 + mu3 * x6])];
        gsym = [0; 0; 0; M \ [1; 0; 0]];
        hsym = (l1 * sin(x1) + l2 * sin(x2) + l3 * sin(x3));

        [f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, degree);

    case 2
        %% Port-Hamiltonian approach to deriving equations of motion
        % redefine Hamiltonian state
        clear x4 x5 x6
        syms p1 p2 p3
        y = [q; p1; p2; p3];

        M = [28/3, 6 * cos(x1 - x2), 2 * cos(x1 - x3);
             6 * cos(x1 - x2), 16/3, 2 * cos(x2 - x3);
             2 * cos(x1 - x3), 2 * cos(x2 - x3), 4/3];
        %         simplify(p - M*qdot) % Check that M satisfies p = M*qdot
        Minv = inv(M); % yikes lol

        T = 1/2 * [p1; p2; p3].' * Minv * [p1; p2; p3];
        H = T + V;

        D = [mu1 + mu2, -mu2, 0;
             -mu2, mu2 + mu3, -mu3;
             0, -mu3, mu3];

        J = [zeros(n / 2), eye(n / 2);
             -eye(n / 2), zeros(n / 2)];
        R = [zeros(n / 2), zeros(n / 2);
             zeros(n / 2), D];

        fsym = (J - R) * gradient(H, y);
        gsym = [0; 0; 0; 1; 0; 0];
        %                 hsym = (l1 * sin(x1) + l2 * sin(x2) + l3/2 * sin(x3));
        hsym = (l1 * sin(x1) + l2 * sin(x2) + l3 * sin(x3));

        [f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, y, degree);
end

end
