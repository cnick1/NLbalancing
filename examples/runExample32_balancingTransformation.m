function runExample32_balancingTransformation(degree,lim,reduction)
%runExample32_balancingTransformation Runs the 2D pendulum example to visualize the nonlinear balancing transformations.
%
%   Usage:  runExample32_balancingTransformation(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description: This simple 3D examples captures the key idea in the
%   model reduction problem: the presence of a subsystem that in some sense
%   contributes little (perhaps is decoupled) to the overall dynamics, yet
%   drives interactions that cannot directly be eliminated.
%       Consider the linear model from [1,2]:
%           ẋ₁ = −x₁ + 100 x₃ + u,
%           ẋ₂ = −2 x₂ + 100 x₃ + u,
%           ẋ₃ = −5 x₃ + u,
%            y = x₁ + x₂ + x₃,
%   The third state component is decoupled and decays quickly, so we
%   intuitively expect that we should be able to approximate this model
%   with a 2D model. However, x₃ strongly drives the states x₁ and x₂. This
%   in some sense directly demonstrates the need for balancing: the state
%   contributes little to observability (since it decays quickly and
%   contributes little to the output) but contributes significantly to
%   controllability (since it drives the other states).
%
%   We compute the energy functions, the input-normal/output-diagonal
%   transformation, and then the true balancing transformation, given by the
%   composition x = ̅Φ(z̄(z̄) = Φ(𝝋(z̄)). The mapping from the z̄ coordinates
%   to the x coordinates is visualized by forming a grid in the z̄ coordinates
%   and mapping that grid to the x coordinates.
%
%   References: [3] P. Holmes, J. L. Lumley, G. Berkooz, and C. W. Rowley,
%                   Turbulence, coherent structures, dynamical systems and
%                   symmetry. Cambridge University Press, 2012. doi:
%                   10.1017/cbo9780511919701.
%
%   Part of the NLbalancing repository.
%%
% close all;
arguments
    degree = 4
    lim = 1
    reduction = true
end
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 32\n')


%% Get system dynamics
[f, g, h] = getSystem32(transform=true);

%% Compute balanced realization
[fbal,gbal,hbal,Tbal] = getBalancedRealization(f,g,h,eta=0,degree=degree-1);
TbalInv = transformationInverse(Tbal);

%% Compute input-normal/output-diagonal realization
[v] = approxPastEnergy(f, g, h, 0, degree);
[w] = approxFutureEnergy(f, g, h, 0, degree);
[~, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree-1);
[finOd,ginOd,hinOd] = transformDynamics(f,g,h,TinOd,degree=degree-1);
[vbal, wbal] = transformEnergyFunctions(v,w,Tbal);
[vinOd, winOd] = transformEnergyFunctions(v,w,TinOd);

fprintf("\n  - FOM dynamics:\n\n")
dispKronPoly(f,degree=degree-1)

fprintf("\n  - Balanced dynamics:\n\n")
dispKronPoly(fbal,degree=degree-1)

fprintf("\n  - Energy Functions:\n\n")
dispKronPoly(v,n=3),dispKronPoly(w,n=3)

fprintf("\n  - Input-normal/output-diagonal energy Functions:\n\n")
dispKronPoly(vinOd,n=3),dispKronPoly(winOd,n=3)

fprintf("\n  - Balanced energy Functions:\n\n")
dispKronPoly(vbal,n=3),dispKronPoly(wbal,n=3)


%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
Fbal = @(z) kronPolyEval(fbal, z);
x0 = [1 1 1].'*.05*lim;

%% Simulate original and transformed systems; same input should give same output
% Solve for z0 initial condition with a Newton type iteration
z0 = kronPolyEval(TbalInv,x0);
fprintf(['\n         -> Initial condition: z0 = [', repmat('%2.2e ', 1, numel(z0)), '], '], z0)
fprintf('       error: %2.2e \n', norm(kronPolyEval(Tbal,z0)-x0))

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
if reduction
    z0 = [z0(1:2); 0];
    Fbaltemp = @(z) kronPolyEval(fbal, z);
    Ir = eye(3); Ir(end) = 0;
    Fbal = @(z) Ir*Fbaltemp(z);
end

% Simulate both systems
[t1, X1] = ode45(@(t, x) F(x), [0, 5], x0);
[t2, Z2] = ode45(@(t, z) Fbal(z), [0, 5], z0);

% Convert z trajectories to x coordinates
X2 = zeros(size(Z2));
for i=1:length(t2)
    X2(i,:) = kronPolyEval(Tbal,Z2(i,:).');
end

% Compute the output for the two trajectories
y1 = zeros(size(t1));
for i=1:length(t1)
    y1(i) = kronPolyEval(h,X1(i,:).');
end
y2 = zeros(size(t2));
for i=1:length(t2)
    y2(i) = kronPolyEval(hbal,Z2(i,:).');
end

% Plot state trajectories
figure('Position', [827 220 560 420]);
subplot(2,1,1); hold on;
plot(t1,X1)
plot(t2,X2,'--')

title('Reconstructed state trajectories'); xlabel('Time t')
ylabel('x_i(t)')
legend('FOM x_1','FOM x_2','FOM x_3','ROM x_1','ROM x_2','ROM x_3')

% Plot outputs
subplot(2,1,2); hold on;
plot(t1,y1)
plot(t2,y2,'--')

title('Model output'); xlabel('Time t')
ylabel('y(t)')
legend('FOM output','ROM output')

fprintf('The output error is: %f \n', norm(interp1(t2, y2, 0:.1:5) - interp1(t1, y1, 0:.1:5)))

%% Manifold figure plot 
dx = 0.05; lim = 2;
[z1,z2,z3] = meshgrid(-lim:dx:lim,-lim:dx:lim,0);
x1 = zeros(size(z1)); x2 = x1; x3 = x1;
for i=1:numel(z1) 
    [x1(i), x2(i), x3(i)] = kronPolyEval(Tbal,[z1(i);z2(i);z3(i)]);
end

figure; 
surf(z1,z2,z3)
figure; 
surf(x1,x2,x3)
drawnow 
return
%% 
global T0; opts = odeset(OutputFcn=@odeprog); T0 = tic;
[T, Z] = ode45(@(t, z) (Az*z + Bz*randn(1)), [0, 50], [0;0;0],opts);
X = zeros(size(Z)); Z = 100*Z;
for i=1:length(T)
    X(i,:) = [Z(i,1), Z(i,2), Z(i,3) - a*(Z(i,1)^2 + Z(i,2)^2) - a*Z(i,1)^3];
end

hold on;
plot3(X(:,1),X(:,2),X(:,3),'r')
xlabel('x_1'); ylabel('x_2'); zlabel('x_3')
title('Balanced manifold and white noise response')

function status = odeprog(t, y, flag)
