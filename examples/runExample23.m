function runExample23(degree,linear)
%runExample23 Runs the 3D example from [1,2].
%
%   Usage:  runExample23(degree)
%
%   Inputs:
%       degree - desired degree of the energy function approximation
%       linear - boolean, whether to use the linear model from [3] or 
%                nonlinear model from [1] (debugging option)\
%
%   Description: These simple 3D examples capture the key idea in the
%   model reduction problem: the presence of a subsystem that in some sense
%   contributes little (perhaps is decoupled) to the overall dynamics, yet
%   drives interactions that cannot directly be eliminated. Consider the
%   linear system from [3]:
%           xdot_1 = −x1 + 100 x3 + u,
%           xdot_2 = −2 x2 + 100 x3 + u,
%           xdot_3 = −5 x3 + u,
%                y = x1 + x2 + x3,
%   The third state component is decoupled and decays quickly, so we
%   intuitively expect that we should be able to apprroximate this model
%   with a 2D model. However, x3 strongly drives the states x1 and x2. This
%   in some sense directly demonstrates the need for balancing: the state
%   contributes little to observability (since it decays quickly and
%   contributes little to the output) but contributes significantly to
%   controllability (since it drives the other states). 
%       The model from [1,2] is a nonlinear example exhibiting the same
%   features. The linear interaction is now replaced with a nonlinear one:
%           xdot_1 = −x1 + 20 x1 x3 + u,
%           xdot_2 = −2 x2 + 20 x2 x3 + u,
%           xdot_3 = −5 x3 + u,
%                y = x1 + x2 + x3,
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
set(groot,'defaultLineLineWidth',1.5)

fprintf('Running Example 23\n')

if nargin < 2
    if nargin < 1
        degree = 10;
    end
    linear = false;
end

[f, g, h] = getSystem23(linear);

%%  Compute the energy functions
fprintf("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~ Computing energy functions:  ~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

eta = 0;
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

% Print some results
n = 3;
x = sym('x', [1, n]).'; syms(x);

fprintf("\n    Controllability energy: \n        Lc = 1/2 *(")
disp(vpa(kronPolyEval(v, x), 3))
fprintf("    Observability energy: \n        Lo = 1/2 *(")
disp(vpa(kronPolyEval(w, x), 3))

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
tic
[sigmaSquared, Tod] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
fprintf("Input-normal/output-diagonal transformation took %f seconds. \n\n", toc)

%% Plot the squared singular value functions
fprintf("    Singular value functions (squared): \n\n")

% figure; hold on;
syms z
for i=1:n
    fprintf("         𝜎_%i^2(z) = tau_%i(z,0) = ",i,i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 4))
    fplot(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 4),[-1 1])
end

set(gca, 'YScale', 'log')

return
%% Compute transformed dynamics
[ft, gt, ht] = transformDynamics(f, g, h, Tod);

% Compare original dynamics with transformed dynamics 
z = sym('z', [1, n]).'; syms(z); syms u;
F = @(x) kronPolyEval(f, x); G = @(x) (g{1} + kronPolyEval(g(2:end), x));
Ft = @(x) kronPolyEval(ft, x); Gt = @(x) (gt{1} + kronPolyEval(gt(2:end), x));

fprintf("\n    Original dynamics: xdot = \n")
disp(vpa(F(x) + G(x)*u, 3))
fprintf("\n    Transformed dynamics: zdot = \n")
disp(vpa(Ft(z) + Gt(z)*u, 3))


%% Simulate original and transformed systems; same input should give same output
x0 = [1 1 1].'/30; 

% [t1, X1] = ode45(@(t, x) F(x), [0, 5], x0*.5);
[t2, X2] = ode45(@(t, x) F(x), [0, 5], x0);

% figure
% subplot(2,1,1)
% plot(t1,X1)
% subplot(2,1,2)
% plot(t1,sum(X1,2))
% 
% figure
% subplot(2,1,1)
% plot(t2,X2)
% subplot(2,1,2)
% plot(t2,sum(X2,2))

figure
hold on
subplot(2,1,1); hold on;
plot(t2,X2)

% plot(t1,sum(X1,2))
subplot(2,1,2); hold on;
plot(t2,sum(X2,2))

z0 = Tod{1}\x0; z0.';
for i=1:10
    z0 = Tod{1}\(x0-kronPolyEval(Tod, z0)+Tod{1}*z0);
    z0.';
    kronPolyEval(Tod, z0);
end
kronPolyEval(Tod, z0)


% [t12, Z1] = ode45(@(t, x) Ft(x), [0, 5], z0*.5);
[t22, Z2] = ode45(@(t, x) Ft(x), [0, 5], z0);

X22 = zeros(size(Z2));
for i=1:length(t22)
    X22(i,:) = kronPolyEval(Tod, Z2(i,:).');
end

subplot(2,1,1)
plot(t22,X22,'--')
subplot(2,1,2)
plot(t22,sum(X22,2),'--')

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
