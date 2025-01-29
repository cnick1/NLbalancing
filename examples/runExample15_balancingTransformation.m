function runExample15_balancingTransformation(degree,reduction,scaling)
%runExample15_balancingTransformation Runs the 4D double pendulum example to visualize the nonlinear balancing transformations.
%
%   Usage:  runExample15_balancingTransformation(degree)
%
%   Inputs:
%       degree    - desired degree of the energy function approximation
%       reduction - boolean, whether or not to apply reduction
%       scaling - factor applied to the initial condition to move it closer
%                   or farther from the origin
%
%   Description: This model has been used several times in the literature [1,2].
%   The 4D state-space model is
%       ẋ = [                  x₃;
%                              x₄;
%              M⁻¹( ∂L/∂(x₁x₂) - Ṁ [x₃; x₄] + Q) ]
%       y = [l₁ * sin(x₁) + l₂ * sin(x₁ + x₂);
%            l₁ * cos(x₁) + l₂ * cos(x₁ + x₂)]
%   The output is taken as the horizontal and vertical positions of the
%   second mass. The mass matrix and its inverse are
%       M(x) = [m₁₁, m₁₂;    M⁻¹(x) = _______1_______  [m₂₂, -m₂₁;
%               m₂₁, m₂₂]            (m₁₁m₂₂ - m₁₂m₂₁) -m₁₂,  m₁₁]
%    where the entries are
%       m₁₁       = m₁ l₁² + m₂ l₁² + m₂ l₂² + 2 m₂ l₁ l₂ cos x₂
%       m₁₂ = m₂₁ = m₂ l₂² + m₂ l₁ l₂ cos x₂
%       m₂₂       = m₂ l₂²
%   The Lagrangian is L(x,ẋ) = T(ẋ) - V(x), where the potential energy is
%       V(x) = - m₁ g l₁ cos x₁ - m₂ g (l₁ cos x₁ + l₂ cos(x₁ + x₂))
%   and the kinetic energy is T(ẋ) = 1/2 ẋ.' * M * ẋ.
%
%   We compute the energy functions, the input-normal/output-diagonal
%   transformation, and then the true balancing transformation, given by the
%   composition x = Φbar(z̄) = Φ(𝝋(z̄)). We then simulate the transformed and
%   optionally reduced system and compare with the original (full-order) model.
%
%   References: [1]  K. Fujimoto and D. Tsubakino, “On computation of
%               nonlinear balanced realization and model reduction,” in
%               2006 American Control Conference, IEEE, 2006. doi:
%               10.1109/acc.2006.1655399
%               [2] K. Fujimoto and D. Tsubakino, “Computation of nonlinear
%               balanced realization and model reduction based on Taylor
%               series expansion,” Systems & Control Letters, vol. 57, no.
%               4, pp. 283–289, Apr. 2008, doi: 10.1016/j.sysconle.2007.08.015
%
%   Part of the NLbalancing repository.
%%
set(groot,'defaultLineLineWidth',1.5,'defaultTextInterpreter','LaTeX')
% clc;
% close all;

fprintf('Running Example 15\n')
if nargin < 3
    scaling = 1;
    if nargin < 2
        reduction = true;
        if nargin < 1
            degree = 6;
        end
    end
end

%% Get system dynamics
[f, g, h] = getSystem15(degree - 1); [p,n] = size(h{1});

%%  Compute the energy functions
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~ Computing energy functions:  ~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
eta = 0;
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
fprintf(" ~~~~~~~~~~~ Computing transformation and singular value functions:  ~~~~~~~~~~~~ \n")
[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
Ft = @(z) PhiBarJacobian(z,TinOd,sigmaSquared)\kronPolyEval(f, PhiBar(z,TinOd,sigmaSquared));

%% Simulate original and transformed systems; same input should give same output
x0 = [0 0 1 -2].'*scaling;

z0 = newtonIteration(x0, TinOd, sigmaSquared);


%% Apply reduction by eliminating z3 (set it and its derivative to zero)
if reduction
    z0(4) = 0; z0(3) = 0;
    Fttemp = @(z) PhiBarJacobian(z,TinOd,sigmaSquared)\kronPolyEval(f, PhiBar(z,TinOd,sigmaSquared));
    Ir = eye(4); Ir(end) = 0; Ir(3,3) = 0;
    Ft = @(z) Ir*Fttemp(z);
end

% Simulate both systems
[t1, X1] = ode45(@(t, x) F(x), [0, 50], x0);
[t2, Z2] = ode45(@(t, z) Ft(z), [0, 50], z0);

% Convert z trajectories to x coordinates
X2 = zeros(size(Z2));
for i=1:length(t2)
    X2(i,:) = PhiBar(Z2(i,:).',TinOd,sigmaSquared);
end

% Get outputs
y1 = zeros(length(t1),p);
for i=1:length(t1)
    y1(i,:) = kronPolyEval(h, X1(i,:).');
end
y2 = zeros(length(t2),p);
for i=1:length(t2)
    y2(i,:) = kronPolyEval(h, X2(i,:).');
end

% Plot state trajectories
figure; hold on;
plot(t1,X1)
plot(t2,X2,'--')

title('Reconstructed state trajectories'); xlabel('Time $t$')
ylim(max(x0)*[-1.5 1.5]); ylabel('$x_i(t)$')
legend('FOM x_1','FOM x_2','FOM x_3','FOM x_4','ROM x_1','ROM x_2','ROM x_3','ROM x_4')

figure; hold on;
plot(t1,y1)
plot(t2,y2,'--')

title('Model output'); xlabel('Time $t$')
ylim(max(max(y1))*[-1.1 1.1]); ylabel('$y(t)$')
legend('FOM output','','ROM output','')


end


