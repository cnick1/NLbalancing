function [f, g, h] = getSystem14(degree, model)
%getSystem14  Generates a polynomial approximation to the 4D model of a double pendulum.
%
%   Usage:  [f,g,h] = getSystem14(degree, relativeAngles)
%
%   References: [1] Tsubakino, [2] Scherpen thesis
%
%%

if nargin < 2
    model = 1;
    if nargin < 1
        degree = 3;
    end
end

n = 4;
x = sym('x', [1, n]).';
syms(x);

G = 9.8;
m1 = 1; m2 = 1;
l1 = 1; l2 = 1;
mu1 = 1; mu2 = 1;

switch model
    case 1 % define joint 2 angle relative to joint 1 angle, as in [1,2]
        P =- m1 * G * l1 * cos(x1) ... % m1 potential energy
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
        Mdot =- [2 * m2 * l1 * l2 * sin(x2) * x4, m2 * l1 * l2 * sin(x2) * x4;
                 m2 * l1 * l2 * sin(x2) * x4, 0];

        K = 0.5 * [x3; x4].' * M * [x3; x4];

        %         D = diag([mu1, mu2]);

        fsym = [x3;
                x4;
                Minv * (gradient(K - P, [x1; x2]) - Mdot * [x3; x4] - [mu1 * x3; mu2 * x4])];
        %             Minv * (gradient(K - P,[x1;x2]) - Mdot * [x3; x4] -  [(mu1+mu2)*x3-mu2*x4;- mu2*x3 + mu2*x4])];

        gsym = [0; 0; Minv * [1; 0]];
        hsym = [l1 * sin(x1) + l2 * sin(x1 + x2);
                l1 * cos(x1) + l2 * cos(x1 + x2)];

        [f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, degree);

        %         g = g(1:3);
        %         h = h(1:3);
    case 2
        % Fujimoto/Scherpen values
        %         m1 = 1; m2 = 10;
        %         l1 = 10; l2 = 1;
        %         G = 9.81;

        % Scherpen/Gray values
        m1 = 1; m2 = 1;
        l1 = 1; l2 = 1;
        G = 10;

        n = 2;
        x = sym('x', [1, n]).';
        syms(x);

        V =- m1 * G * l1 * cos(x1) ... % m1 potential energy
            - m2 * G * (l1 * cos(x1) + l2 * cos(x1 + x2)); % m2 potential energy

        m11 = m1 * l1 ^ 2 + m2 * l1 ^ 2 + m2 * l2 ^ 2 + 2 * m2 * l1 * l2 * cos(x2);
        m12 = m2 * l2 ^ 2 + m2 * l1 * l2 * cos(x2);
        m22 = m2 * l2 ^ 2;

        M = [m11, m12;
             m12, m22];
        Minv = 1 / (m11 * m22 - m12 ^ 2) * ...
            [m22, -m12;
         -m12, m11];

        fsym = [-Minv * gradient(V, x)];
        gsym = [Minv * [1; 0]];
        hsym = x1;

        [f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, degree);

        %         g = full(g{1}); h = full(h{1});

    case 3
        % Fujimoto/Scherpen values
        m1 = 1; m2 = 10;
        l1 = 10; l2 = 1;
        G = 9.81;

        n = 2;
        x = sym('x', [1, n]).';
        syms(x);

        V =- m1 * G * l1 * cos(x1) ... % m1 potential energy
            - m2 * G * (l1 * cos(x1) + l2 * cos(x1 + x2)); % m2 potential energy

        m11 = m1 * l1 ^ 2 + m2 * l1 ^ 2 + m2 * l2 ^ 2 + 2 * m2 * l1 * l2 * cos(x2);
        m12 = m2 * l2 ^ 2 + m2 * l1 * l2 * cos(x2);
        m22 = m2 * l2 ^ 2;

        M = [m11, m12;
             m12, m22];
        Minv = 1 / (m11 * m22 - m12 ^ 2) * ...
            [m22, -m12;
         -m12, m11];

        fsym = [-Minv * gradient(V, x)];
        gsym = [Minv * [1; 0]];
        hsym = x1;

        [f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, degree);

        %         g = full(g{1}); h = full(h{1});

    case 4
        % define positions from fixed vertical, see e.g. https://scienceworld.wolfram.com/physics/DoublePendulum.html
        % Using velocities as x3 and x4
        % T = 0.5 * m1 * l1 ^ 2 * x3 ^ 2 ...                                                           % m1 kinetic energy
        %   + 0.5 * m2 * (l1 ^ 2 * x3 ^ 2 + l2 ^ 2 * x4 ^ 2 + 2 * l1 * l2 * x3 * x4 * cos(x1 - x2));   % m2 kinetic energy

        % Using momenta as x3 and x4
        Minv = 1 / (m2 * l1 ^ 2 * l2 ^ 2 * (m1 + m2 * sin(x1 - x2) ^ 2)) ...
            * [m2 * l2 ^ 2, -m2 * l1 * l2 * cos(x1 - x2);
           -m2 * l1 * l2 * cos(x1 - x2), l1 ^ 2 * (m1 + m2)];

        T = 0.5 * [x3; x4].' * Minv * [x3; x4];
        V =- m1 * G * l1 * cos(x1) ... % m1 potential energy
            - m2 * G * (l1 * cos(x1) + l2 * cos(x2)); % m2 potential energy

        Hamiltonian = T + V;

        D = diag([mu1, mu2]);

        J = [zeros(n / 2), eye(n / 2);
             -eye(n / 2), zeros(n / 2)];
        R = [zeros(n / 2), zeros(n / 2);
             zeros(n / 2), D];

        fsym = (J - R) * gradient(Hamiltonian, x);
        gsym = [0; 0; 1; 0];
        hsym = [l1 * sin(x1) + l2 * sin(x1 + x2);
                l1 * cos(x1) + l2 * cos(x1 + x2)];

        [f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, degree);

        g = full(g{1}); h = full(h{1});
end

end
