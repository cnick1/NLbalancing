function [T1, T2, T3] = runExample3_balancedReduction(degree, n, r,  U0, verbose)
%runExample3_balancedReduction 
%
%   Usage:  runExample3_balancedReduction(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description: 
%
%   References:
%
%   Part of the NLbalancing repository.
%%
% close all;
arguments
    degree = 2
    n = 33
    r = 8
    U0 = 0
    verbose = false
end
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 3, polynomial balanced reduced-order model...\n')

close all

%% Get dynamics and define control problem
[f, g, h, initialCondition] = getSystem3(n, 4, 4, 0.001, 0);
xg = linspace(0,1,n);

g = {diag(-4*(xg-.5).^2+1)};
h = {eye(n)};

x0 = 1e-1*initialCondition;

if verbose
    fprintf('  - The full-order nonlinear model is:\n')
    dispPolyDynamics(f,g,h)
end

%% Compute balanced realization
if verbose
    [fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g,h,eta=0,degree=3,transformationDegree=degree-1);
    fprintf('  - The full-order balanced realization for the nonlinear model is:\n')
    dispPolyDynamics(fbal,gbal,hbal,variable='z')
end

tic; fprintf('   Computing balanced reduced-order model...')
[fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g,h,r=r,eta=0,degree=3,transformationDegree=degree-1);
TbalInv = transformationInverse(Tbal);
T1 = toc; fprintf("completed in %2.2f seconds. \n", T1)

if verbose
    plotSingularValueFunctions(sigmaSquared)
    fprintf('  - The reduced-order balanced realization for the nonlinear model is:\n')
    dispPolyDynamics(fbal,gbal,hbal,variable='z')
end

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
% G = @(x) kronPolyEval(g, x, scenario='G(x)');
G = @(x) kronPolyEval(g, x, scenario='G(x)');
Fbal = @(z) kronPolyEval(fbal, z);
% Gbal = @(z) kronPolyEval(gbal, z, scenario='G(x)');
Gbal = @(z) kronPolyEval(gbal, z, scenario='G(x)');
m = size(g{1},2);
u = @(t) [zeros((m-1)/2,1); U0*sin(2*pi*t/5); zeros((m-1)/2,1)];

%% Simulate original and transformed systems; same input should give same output
% Solve for z0 initial condition
z0 = kronPolyEval(TbalInv,x0);
% z0 = newtonIteration(x0, @(z) kronPolyEval(Tbal,z), @(z) jcbn(Tbal,z),maxIter=100,verbose=true,z0=kronPolyEval(TbalInv,x0));

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
% Simulate both systems
global T0; 
opts = odeset(OutputFcn=@odeprog);
tmax = 10; t = 0:.1:tmax; % specify for plotting

fprintf('   Simulating the FOM in the original coordinates...')
T0 = tic; [t1, X1] = ode15s(@(t, x) F(x) + G(x)*u(t), t, x0, opts);
T2 = toc(T0); fprintf("completed in %2.2f seconds. \n", T2)


fprintf('   Simulating the ROM in the balanced coordinates...')
T0 = tic; [t2, Z2] = ode15s(@(t, z) Fbal(z) + Gbal(z)*u(t), t, z0, opts);
T3 = toc(T0); fprintf("completed in %2.2f seconds. \n", T3)

%% Compute and plot the outputs 
p = size(h{1},1);
y1 = zeros(length(t1),p);
y2 = zeros(length(t2),p);
for i=1:length(t1) 
    y1(i,:) = kronPolyEval(h, X1(i,:).');
end
for i=1:length(t2) 
    y2(i,:) = kronPolyEval(hbal, Z2(i,:).');
end

%% Plot solution
figure('Position',[501 19 632 638])
% figure('Position',[474.3333 340.3333 925.3333 300.6667]);
subplot(2,1,1)
mesh(xg,t1,y1);
grid on, axis([0 1 0 tmax -.25e-3 .25e-3]), view(145,35), colormap([0 0 0]);
xlabel x, ylabel t, zlabel u(x,t), title FOM; drawnow

% figure('Position',[474.3333 340.3333 925.3333 300.6667]);
subplot(2,1,2)
mesh(xg,t2,y2);
grid on, axis([0 1 0 tmax -.25e-3 .25e-3]), view(145,35), colormap([0 0 0]);
xlabel x, ylabel t, zlabel u(x,t), title ROM; drawnow

figure
plot(t1,y1(:,p-2))
hold on;
plot(t2,y2(:,p-2))


end

