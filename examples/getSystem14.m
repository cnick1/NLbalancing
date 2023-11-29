function [f, g, h] = getSystem14(degree, model)
%getSystem14  Returns a polynomial approximation to the 2D gradient model of a double pendulum.
%
%   Usage:  [f,g,h] = getSystem14(degree, model)
%
%   Inputs:
%       degree - desired degree of the computed polynomial dynamics.
%       model  - used to select either the Scherpen/Gray [1,2] model (1 - default)
%                or the Fujimoto [3] model (2).
%                In [1,2], the values used are m1=m2=l1=l2=1 and G = 10;
%                in [3], the values used are m1=l2=10, l1=m2=1 and G = 9.81.
%
%   Outputs:
%       f,g,h  - Cell arrays containing the polynomial coefficients for the
%                drift, input, and output (generalizations containing A,B,C)
%
%   Background: This model has been used several times in the literature [1-3].
%       The system describes a set of dynamics related to the double
%       pendulum; however, where the double pendulum would have 4D dynamics
%       and only marginal stability, the associated 2D gradient system is
%       asymptotically stable, and hence the model was more approachable
%       when method were more limited.
%
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
%       The full dynamics are 4 dimensional and are not asymptotically
%       stable. However, the gradient system dynamics are 2D and
%       asymptotically stable. The gradient system dynamics are
%           xdot = -M^{-1}(x) 𝜕V(x)/𝜕x + M^{-1}(x)[1;0] u
%           y    = x1
%
%   References: [1] J. M. A. Scherpen, “Balancing for nonlinear systems,”
%               PhD Dissertation, University of Twente, 1994.
%               [2] W. S. Gray and J. M. A. Scherpen, “Minimality and local
%               state decompositions of a nonlinear state space realization
%               using energy functions,” IEEE Transactions on Automatic
%               Control, vol. 45, no. 11, pp. 2079–2086, 2000, doi:
%               10.1109/9.887630
%               [3] K. Fujimoto and J. M. A. Scherpen, “Nonlinear
%               input-normal realizations based on the differential
%               eigenstructure of Hankel operators,” IEEE Transactions on
%               Automatic Control, vol. 50, no. 1, pp. 2–18, Jan. 2005,
%               doi: 10.1109/tac.2004.840476
%
%
%%

if nargin < 2
    model = 1;
    if nargin < 1
        degree = 3;
    end
end

switch model
    case 1
        % Scherpen/Gray values
        m1 = 1; m2 = 1;
        l1 = 1; l2 = 1;
        G = 10;
    case 2
        % Fujimoto/Scherpen values
        m1 = 1; m2 = 10;
        l1 = 10; l2 = 1;
        G = 9.81;
end

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
Minv = 1 / (m11 * m22 - m12 ^ 2) * [m22, -m12;
                                     -m12, m11];

fsym = -Minv * gradient(V, x);
gsym = Minv * [1; 0];
hsym = x1;

[f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, degree);

end
