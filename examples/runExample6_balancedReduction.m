function [T1, T2, T3, error] = runExample6_balancedReduction(degree, numEls, r,  U0, verbose, nvp)
%runExample6_balancedReduction 
%
%   Usage:  runExample6_balancedReduction(degree,lim)
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
    degree = 4
    numEls = 4
    r = 6*(numEls+1)
    U0 = 1.5e-2
    verbose = false
    nvp.plot = true
end
set(groot,'defaultLineLineWidth',1.5,'defaultTextInterpreter','latex')

numNodes = numEls + 1; n = 6*numEls;
fprintf('Running Example 6, n=%d, degree %d...\n',n,degree)

%% Get dynamics and define control problem
[E, f, g, ~, initialCondition] = getSystem6(numEls, 3, false, false);
g = g(1); %g{1} = g{1}(:,1);
h = {[eye(n/2), zeros(n/2)]};

g = {eye(n)};
h = {eye(n)};

x0 = 0*initialCondition;

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
G = @(x) kronPolyEval(g, x, scenario='G(x)');
Fbal = @(z) kronPolyEval(fbal, z);
Gbal = @(z) kronPolyEval(gbal, z, scenario='G(x)');
m = size(g{1},2);
u = @(t) [zeros(m-1,1); U0*sin(2*pi*t/10)];

%% Simulate original and transformed systems; same input should give same output
% Solve for z0 initial condition
z0 = kronPolyEval(TbalInv,x0);
% z0 = newtonIteration(x0, @(z) kronPolyEval(Tbal,z), @(z) jcbn(Tbal,z),maxIter=100,verbose=true,z0=kronPolyEval(TbalInv,x0));

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
% Simulate both systems
global T0; 
opts = odeset(OutputFcn=@odeprog);

fprintf('   Simulating the FOM in the original coordinates...')
T0 = tic; [t1, X1] = ode23s(@(t, x) F(x) + G(x)*u(t), [0:.1:5], x0, opts);
T2 = toc(T0); fprintf("completed in %2.2f seconds. \n", T2)


fprintf('   Simulating the ROM in the balanced coordinates...')
T0 = tic; [t2, Z2] = ode23s(@(t, z) Fbal(z) + Gbal(z)*u(t), [0:.1:5], z0, opts);
T3 = toc(T0); fprintf("completed in %2.2f seconds. \n", T3)

%% Compute and plot the outputs 
p = size(h{1},1);
y1 = zeros(p,length(t1));
y2 = zeros(p,length(t2));
for i=1:length(t1) 
    y1(:,i) = kronPolyEval(h, X1(i,:).');
end
for i=1:length(t2) 
    y2(:,i) = kronPolyEval(hbal, Z2(i,:).');
end

% figure
% plot(y1(n/2-2,:))
% hold on
% plot(y2(n/2-2,:))

if n == 6 || nvp.plot
    figure(21414836)
    if degree == 2
        close
        figure(21414836)
        plot(t1,y1(n/2-2,:),'DisplayName','FOM output')
        hold on
        plot(t2,y2(n/2-2,:),'--','DisplayName','ROM output w/ linear transformation')
    elseif degree == 3
        plot(t2,y2(n/2-2,:),':o','DisplayName',sprintf('ROM output w/ degree %i transformation',degree-1))
    else
        plot(t2,y2(n/2-2,:),'--+','DisplayName',sprintf('ROM output w/ degree %i transformation',degree-1))
        xlabel('Time, t'); ylabel('Beam tip horizontal displacement $y(t)$');
        legend('Location','southwest'); %ylim(1e-7*[-.35 .05])
        set(gcf,"Position", [573 443.6667 720 400])
        exportgraphics(gca, sprintf('plots/example6_n%i_d%i_U0%i_y.pdf',n,degree,U0*100), 'ContentType', 'vector');
    end
end

error = norm(interp1(t2, y2(n/2-2,:), 0:.1:5) - interp1(t1, y1(n/2-2,:), 0:.1:5));
fprintf('\n   The output error is: ||yᵣ(t)-y(t)||₂ = %e \n\n', error)

end

