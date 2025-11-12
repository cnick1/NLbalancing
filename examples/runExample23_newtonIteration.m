function runExample23_newtonIteration(degree,scaling,reduction)
%runExample23_newtonIteration Runs the Otto 3D example to compare linear vs nonlinear balancing.
%
%   Usage:  runExample23_newtonIteration(degree,reduction,scaling)
%
%   Inputs:
%       degree    - desired degree of the energy function approximation
%       scaling   - factor applied to the initial condition to move it closer
%                   or farther from the origin
%       reduction - boolean, whether or not to apply reduction
%
%   Description: This simple 3D example captures the key idea in the
%   model reduction problem: the presence of a subsystem that in some sense
%   contributes little (perhaps is decoupled) to the overall dynamics, yet
%   drives interactions that cannot directly be eliminated.
%       Consider the nonlinear model from [1,2]:
%           ẋ₁ = −x₁ + 20 x₁ x₃ + u,
%           ẋ₂ = −2 x₂ + 20 x₂ x₃ + u,
%           ẋ₃ = −5 x₃ + u,
%            y = x₁ + x₂ + x₃,
%   The third state component is decoupled and decays quickly, so we
%   intuitively expect that we should be able to approximate this model
%   with a 2D model. However, x₃ strongly drives the states x₁ and x₂. This
%   in some sense directly demonstrates the need for balancing: the state
%   contributes little to observability (since it decays quickly and
%   contributes little to the output) but contributes significantly to
%   controllability (since it drives the other states). This model is a
%   nonlinear adaptation of the linear balancing example Holmes/Rowley
%   showed in [3], see runExample23_holmes().
%
%   Despite being so simple, this is a challenging problem because the
%   nonlinear interaction is strong: those terms are much larger than the
%   linear terms!
%
%   In this script, we use linear vs nonlinear balancing to try to approximate
%   the same test cases from [1-3], i.e. approximating the impulse response
%   of these systems with a 2D ROM.
%
%   Note: We can also look at just the balancing transformation without
%   reduction. In this case, if the transformation is valid (bijective) and
%   well conditioned, the transformed system should be identical to the
%   original system. However, we know that the nonlinear transformations
%   are only valid transformations locally; on the other hand, the linear
%   transformation may only balance the system locally, but it is a valid
%   transformation globally!
%
%   The recommended test cases are:
%       runExample23(2)
%       runExample23(4) % fails due to singularities
%
%       runExample23(2,0.05)
%
%   References: [1] S. E. Otto, A. Padovan, and C. W. Rowley, "Optimizing
%                   oblique projections for nonlinear systems using
%                   trajectories," SIAM Journal on Scientific Computing,
%                   vol. 44, no. 3, pp. A1681–A1702, Jun. 2022, doi:
%                   10.1137/21m1425815.
%               [2] S. E. Otto, A. Padovan, and C. W. Rowley, "Model
%                   reduction for nonlinear systems by balanced truncation
%                   of state and gradient covariance,” SIAM Journal on
%                   Scientific Computing, vol. 45, no. 5, pp. A2325–A2355,
%                   Sep. 2023, doi: 10.1137/22m1513228.
%               [3] P. Holmes, J. L. Lumley, G. Berkooz, and C. W. Rowley,
%                   Turbulence, coherent structures, dynamical systems and
%                   symmetry. Cambridge University Press, 2012. doi:
%                   10.1017/cbo9780511919701.
%
%   Part of the NLbalancing repository.
%%
set(groot,'defaultLineLineWidth',1.5,'defaultTextInterpreter','LaTeX')
% clc;
% close all;

fprintf('Running Example 23\n')

if nargin < 3
    reduction = true;
    if nargin < 2
        scaling = 0.5;
        if nargin < 1
            degree = 2;
        end
    end
end

%% Get system dynamics
[f, g, h] = getSystem23();
f{2} = f{2}*1/10;

%%  Compute the energy functions
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~ Computing energy functions:  ~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
eta = 0;
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
fprintf(" ~~~~~~~~~~~ Computing transformation and singular value functions:  ~~~~~~~~~~~~ \n")
[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree - 1);

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
Ft = @(z) PhiBarJacobian(z,TinOd,sigmaSquared)\kronPolyEval(f, PhiBar(z,TinOd,sigmaSquared));

%% Simulate original and transformed systems; same input should give same output
x0 = [1 1 1].'*scaling;

% Solve for z0 initial condition with a Newton type iteration
z0 = newtonIteration(x0, @(z) PhiBar(z,TinOd,sigmaSquared), @(z) PhiBarJacobian(z,TinOd,sigmaSquared),maxIter=100,verbose=true);

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
if reduction
    z0(3) = 0;
    Fttemp = @(z) PhiBarJacobian(z,TinOd,sigmaSquared)\kronPolyEval(f, PhiBar(z,TinOd,sigmaSquared));
    Ir = eye(3); Ir(end) = 0;
    Ft = @(z) Ir*Fttemp(z);
end

% Simulate both systems
[t1, X1] = ode45(@(t, x) F(x), [0, 5], x0);
[t2, Z2] = ode45(@(t, z) Ft(z), [0, 5], z0);

% Convert z trajectories to x coordinates
X2 = zeros(size(Z2));
for i=1:length(t2)
    X2(i,:) = PhiBar(Z2(i,:).',TinOd,sigmaSquared);
end

% Plot state trajectories
figure('Position', [827 220 560 420]);
subplot(2,1,1); hold on;
plot(t1,X1)
plot(t2,X2,'--')

title('Reconstructed state trajectories'); xlabel('Time $t$')
ylabel('$x_i(t)$')
legend('FOM $x_1$','FOM $x_2$','FOM $x_3$','ROM $x_1$','ROM $x_2$','ROM $x_3$')

% Plot outputs
subplot(2,1,2); hold on;
plot(t1,sum(X1,2))
plot(t2,sum(X2,2),'--')

title('Model output'); xlabel('Time $t$')
ylabel('$y(t)$')
legend('FOM output','ROM output')

fprintf('The output error is: %f \n', norm(interp1(t2, sum(X2,2), 0:.1:5) - interp1(t1, sum(X1,2), 0:.1:5)))

end


