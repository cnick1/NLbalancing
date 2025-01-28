function runExample15_balancingTransformation(degree)
%runExample15_balancingTransformation Runs the 4D double pendulum example to visualize the nonlinear balancing transformations.
%
%   Usage:  runExample15_balancingTransformation(degree)
%
%   Inputs:
%       degree    - desired degree of the energy function approximation
%
%   Description: 
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
scaling = 1; reduction = true;
if nargin < 1
    degree = 6;
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


