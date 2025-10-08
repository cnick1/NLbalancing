function [f, g, h] = getSystem26(degree,transformed)
%getSystem26 Generates a polynomial 2D inverted pendulum model
%
%   Usage:  [f,g,h] = getSystem26()
%
%   Inputs:    degree - desired degree of the polynomial approximation
%         transformed - whether to use model (1) or (2)
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The basic equation for the inverted pendulum is
%
%             ẍ - sin(x) = u(t)
%
%   This can be put in first-order form with x₁ = x, x₂ = ẋ as
%
%             ẋ₁ = x₂
%      (1)    ẋ₂ = sin(x₁) + u(t)
%              y = x₁
%
%   Alternatively, we can transform with the following (bijective)
%   coordinate transformation
%
%              x = Φ(z)   = [z₁; z₂ - 2 sin(z₁/2)]
%              z = Φ⁻¹(x) = [x₁; x₂ + 2 sin(x₁/2)]
%
%   to obtain the transformed dynamical system
%
%             ̇z₁ = z₂ - 2 sin(z₁/2)
%      (2)    ̇z₂ = z₂ cos(z₁/2) + u(t)
%              y = z₁
%
%   ( This is done by taking the original dynamics ẋ = f(x) + g(x) u and
%     transforming to get ∂Φ(z)/∂z ż = f(Φ(z)) + g(Φ(z)) u, which upon
%     inverting the Jacobian gives ż =[∂Φ(z)/∂z]⁻¹ f(Φ(z)) + g(Φ(z)) u . )
%
%   More background: The pendulum is a Hamiltonian system. Open-loop
%   trajectories correspond to level-sets of the Hamiltonian, since energy
%   is conserved. Defining q = θ and ̇q = dθ/dt, the kinetic energy is T =
%   ½ m r² ̇q² , whereas the potential energy is V = m g r (cos(q) - 1). If
%   we set all the constants to 1, the Lagrangian is
%
%           L(q, ̇q) = T - V = ½ ̇q² - (cos(q)-1)
%                           = ½ ̇q² - cos(q) + 1
%
%   The Euler-Lagrange equation is d/dt(𝜕L/𝜕q̇) - 𝜕L/𝜕q  = 0, which yields
%   the second-order form of the governing equation
%
%             q̈ - sin(q) = 0
%
%   The Hamiltonian instead is the total energy H = T + V; level-sets of
%   the Hamiltonian are curves of constant energy, so open-loop
%   trajectories which conserve energy are level-sets of the Hamiltonian.
%   So we can compute the separatrix as a particular level set; the way we
%   set up the potential and kinetic energy, the equilibrium is at T = 0
%   and V = 0, so the separatrix is the curve H = 0, which is given by
%
%           ½ ̇q² + (cos(q)-1) = 0    →    ̇q = ± √(2(1-cos(q)))
%
%   These curves are identical to ̇q = ± 2 sin(q/2). Now, in the standard
%   state-space coordinates, we define x₁ = q, x₂ = ̇q, so the curve we are
%   interested in is x₂ = - 2 sin(x₁/2), which we wish to transform to lie
%   along the horizontal axis in the new coordinates, i.e. z₂ = 0. So we
%   want z₂ = 0 = x₂ + 2 sin(x₁/2). This is how we get the transformation
%
%              x = Φ(z)   = [z₁; z₂ - 2 sin(z₁/2)]
%              z = Φ⁻¹(x) = [x₁; x₂ + 2 sin(x₁/2)]
%
%   The Jacobian of this transformation and its inverse are
%
%              [∂Φ(z)/∂z]   = [1 0; -cos(z₁/2)  1]
%              [∂Φ(z)/∂z]⁻¹ = [1 0;  cos(z₁/2)  1]
%
%             ẋ₁ = x₂
%             ẋ₂ = sin(x₁) + u(t)
%
%             ẋ₁ = z₂ - 2 sin(z₁/2)
%             ẋ₂ = sin(z₁) + u(t)
%
%             ̇z₁ = z₂ - 2 sin(z₁/2)
%             ̇z₂ = cos(z₁/2)(z₂ - 2 sin(z₁/2)) + sin(z₁) + u(t)
%                 = z₂cos(z₁/2) - 2 cos(z₁/2) sin(z₁/2) + sin(z₁) + u(t)
%                 = z₂cos(z₁/2) - sin(z₁) + sin(z₁) + u(t)
%                 = z₂cos(z₁/2) + u(t)
%
%
%   References: [1]
%
%
%%
if nargin < 2
    transformed = false;
    if nargin < 1
        degree = 9;
    end
end

n = 2;

if transformed
    %       ̇z₁ = -2 sin(z₁/2) + z₂
    %  (2)  ̇z₂ = z₂ cos(z₁/2) + u(t)
    %        y = z₁
    
    A = [-1 1;
        0 1];
    F2 = sparse(n,n^2);
    f = {A, F2};
    for i = 1:(degree - 1) / 2
        f{end + 1} = sparse(2, 2 ^ (2 * i + 1));
        
        % form ̇z₁ = -2 sin(z₁/2) + z₂
        f{end}(1, 1) = -(-1) ^ i * (1/2)^(2 * i + 1) / factorial(2 * i + 1) * 2; % -2 sin(z₁/2)
        
        % form ̇z₂ = z₂ cos(z₁/2) + u(t)
        f{end}(2, 2) = (-1) ^ i * (1/2)^(2 * i) / factorial(2 * i);          % z₂ cos(z₁/2)
        
        f{end + 1} = sparse(2, 2 ^ (2 * i + 2));
    end
    f = f(1:degree);
    
    % Symbolic approach
    % z = sym('z', [1, 2]).'; syms(z);
    % % % % fsym = [1 0; 1/4*cos(z1/2) 1]*[z2 - 1/2*sin(z1/2); sin(z1)]; % Sanity check, identical
    % fsym = [1 0; -cos(z1/2) 1]\[z2 - 2*sin(z1/2); sin(z1)]; % Sanity check, identical
    % gsym = [0;1]; hsym = z1;
    % [f2,~,~] = approxPolynomialDynamics(fsym,gsym,hsym,z,degree);
    
    
else % standard state-space form
    %        ẋ₁ = x₂
    %   (1)  ẋ₂ = sin(x₁) + u(t)
    %         y = x₁
    
    A = [0 1;
        1 0];
    F2 = sparse(n,n^2);
    f = {A, F2};
    for i = 1:(degree - 1) / 2
        f{end + 1} = sparse(2, 2 ^ (2 * i + 1));
        f{end}(2, 1) = (-1) ^ i / factorial(2 * i + 1);
        f{end + 1} = sparse(2, 2 ^ (2 * i + 2));
    end
    f = f(1:degree);
    
end

B = [0; 1];
C = [1; 0];

g = {B};
h = {C};

end
