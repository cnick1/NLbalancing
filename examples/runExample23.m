function runExample23(degree,linear,reduction,scaling)
%runExample23 Runs the 3D example from [1,2].
%
%   Usage:  runExample23(degree)
%
%   Inputs:
%       degree    - desired degree of the energy function approximation
%       linear    - boolean, whether to use the linear model from [3] or
%                   nonlinear model from [1] (debugging option)
%       reduction - boolean, whether or not to apply reduction
%       scaling - factor applied to the initial condition to move it closer
%                   or farther from the origin
%
%   Description: These simple 3D examples capture the key idea in the
%   model reduction problem: the presence of a subsystem that in some sense
%   contributes little (perhaps is decoupled) to the overall dynamics, yet
%   drives interactions that cannot directly be eliminated. Consider the
%   linear system from [3]:
%           ẋ₁ = −x₁ + 100 x₃ + u,
%           ẋ₂ = −2 x₂ + 100 x₃ + u,
%           ẋ₃ = −5 x₃ + u,
%                y = x₁ + x₂ + x₃,
%   The third state component is decoupled and decays quickly, so we
%   intuitively expect that we should be able to apprroximate this model
%   with a 2D model. However, x₃ strongly drives the states x₁ and x₂. This
%   in some sense directly demonstrates the need for balancing: the state
%   contributes little to observability (since it decays quickly and
%   contributes little to the output) but contributes significantly to
%   controllability (since it drives the other states).
%       The model from [1,2] is a nonlinear example exhibiting the same
%   features. The linear interaction is now replaced with a nonlinear one:
%           ẋ₁ = −x₁ + 20 x₁ x₃ + u,
%           ẋ₂ = −2 x₂ + 20 x₂ x₃ + u,
%           ẋ₃ = −5 x₃ + u,
%                y = x₁ + x₂ + x₃,
%   Despite being so simple, this is a challenging problem because the
%   nonlinear interaction is strong: those terms are much larger than the
%   linear terms!
%
%   In this script, we use nonlinear balancing to try to approximate the
%   same test cases from [1,2,3], i.e. approximating the impulse response
%   of these systems with a 2D ROM.
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

if nargin < 4
    scaling = 1/20;
    if nargin < 3
        reduction = false;
        if nargin < 2
            linear = false;
            if nargin < 1
                if linear
                    degree = 2;
                else
                    degree = 4;
                end
            end
        end
    end
end

[f, g, h] = getSystem23(linear);

%%  Compute the energy functions
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~ Computing energy functions:  ~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

eta = 0;
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

% Print energy functions
n = 3; x = sym('x', [1, n]).'; syms(x);

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
fprintf(" ~~~~~~~~~~~ Computing transformation and singular value functions:  ~~~~~~~~~~~~ \n")
[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree - 1);

%% Plot the squared singular value functions
fprintf("    Singular value functions (squared): \n")
figure('Position', [266 220 560 420]); hold on; syms z
for i=1:n
    fprintf("         𝜎_%i^2(z) = tau_%i(z,0) = ",i,i)
    fprintf(string(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 4)));fprintf("\n")
    fplot(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 4),[-1 1],'LineWidth',1.5)
end

title('Singular value functions (squared)')
xlabel('$z_i$'); ylabel('$\sigma_i^2$'); set(gca, 'YScale', 'log')
legend('\sigma_1^2','\sigma_2^2','\sigma_3^2')

%% Apply reduction by eliminating all appearances of z3 from the transformation???
% TinOd{1} = TinOd{1}(:,1:2);

%% Compute transformed dynamics
[ft, gt, ht] = transformDynamics(f, g, h, TinOd);

% Compare original dynamics with transformed dynamics
z = sym('z', [1, n]).'; syms(z); syms u;
F = @(x) kronPolyEval(f, x); G = @(x) (g{1} + kronPolyEval(g(2:end), x));
Ft = @(z) kronPolyEval(ft, z); Gt = @(z) (gt{1} + kronPolyEval(gt(2:end), z));

%% Simulate original and transformed systems; same input should give same output
x0 = [1 1 1].'*scaling;

% Solve for z0 initial condition with a Newton type iteration
z0 = newtonIteration(x0, @(z) kronPolyEval(TinOd, z), @(z) jcbn(TinOd, z));

% Simulate both systems
[t2, X2] = ode45(@(t, x) F(x), [0, 5], x0);
[t22, Z2] = ode45(@(t, x) Ft(x), [0, 5], z0);

% Convert z trajectories to x coordinates
X22 = zeros(size(Z2));
for i=1:length(t22)
    X22(i,:) = kronPolyEval(TinOd, Z2(i,:).');
end

% Plot state trajectories
figure('Position', [827 220 560 420]);
subplot(2,1,1); hold on;
plot(t2,X2)
plot(t22,X22,'--')

title('Reconstructed state trajectories'); xlabel('Time $t$')
ylim(max(x0)*[-1.5 1.5]); ylabel('$x_i(t)$')
legend('FOM x_1','FOM x_2','FOM x_3','ROM x_1','ROM x_2','ROM x_3')

% Plot outputs
subplot(2,1,2); hold on;
plot(t2,sum(X2,2))
plot(t22,sum(X22,2),'--')

title('Model output'); xlabel('Time $t$')
ylim(max(x0)*[-3.5 3.5]); ylabel('$y(t)$')
legend('FOM output','ROM output')

%% Truncate last state to obtain a reduced-order model
% for i=1:length(ft)
%     fr{i} = ft{i}(1:2,:);
% end
%
% for i=1:length(gt)
%     fr{i} = ft{i}(1:2,:);
% end
%
% Fr = @(x) kronPolyEval(fr, x); Gr = @(x) (gr{1} + kronPolyEval(gr(2:end), x));



end


