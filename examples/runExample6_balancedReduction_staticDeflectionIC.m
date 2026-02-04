function [T1, T2, T3, error, x0] = runExample6_balancedReduction_staticDeflectionIC(x0, degree, numEls, r,  U0, verbose, figNum, nvp)
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
    x0 = []
    degree = 4
    numEls = 1
    r = 4
    U0 = 3e4
    verbose = false
    figNum = 21414836
    nvp.plot = true
end
set(groot,'defaultLineLineWidth',2,'defaultTextInterpreter','latex')

numNodes = numEls + 1; n = 6*numEls;
fprintf('Running Example 6, n=%d, degree %d...\n',n,degree)

%% Get dynamics and define control problem
[E, f, g, ~, initialCondition] = getSystem6(numEls, 3, false, false);

g = {eye(n)}; h = {eye(n)}; m = size(g{1},2); p = size(h{1},1);
F = @(x) kronPolyEval(f, x);
G = @(x) kronPolyEval(g, x, scenario='G(x)');

if isempty(x0)
    % Option 1: Obtain steady-state via time integration
    % u = @(t) [zeros(m-2,1); U0*(t>0.0005); 0];
    % [t, X1] = ode23s(@(t, x) F(x) + G(x)*u(t), [0 .1], 0*initialCondition);
    % y1 = zeros(p,length(t));
    % for i=1:length(t)
    %     y1(:,i) = kronPolyEval(h, X1(i,:).');
    % end
    % figure
    % plot(y1(n/2-2,:))
    % 
    % x0 = X1(end,:).';


    % Option 2: Obtain steady-state Newton iteration for equilibrium point
    fsymmetric = f; 
    for i=2:length(fsymmetric)
        for j = 1:n
            fsymmetric{i}(j,:) = kronMonomialSymmetrize(fsymmetric{i}(j,:),n,i);
        end
    end

    u = [zeros(m-2,1); U0; 0]; 
    x0 = newtonIteration(-g{1}*u, @(x) kronPolyEval(f, x), @(x) jcbn(fsymmetric, x), maxIter=100);
end

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

Fbal = @(z) kronPolyEval(fbal, z);
Gbal = @(z) kronPolyEval(gbal, z, scenario='G(x)');
u = @(t) [zeros(m,1)];

%% Simulate original and transformed systems; same input should give same output
% Solve for z0 initial condition
z0 = kronPolyEval(TbalInv,x0);
% z0 = newtonIteration(x0, @(z) kronPolyEval(Tbal,z), @(z) jcbn(Tbal,z),maxIter=100,verbose=true,z0=kronPolyEval(TbalInv,x0));

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
% Simulate both systems
global T0;
opts = odeset(OutputFcn=@odeprog);
t = [0:.0001:.05];

fprintf('   Simulating the FOM in the original coordinates...')
T0 = tic; [t, X1] = ode23s(@(t, x) F(x) + G(x)*u(t), t, x0, opts);
T2 = toc(T0); fprintf("completed in %2.2f seconds. \n", T2)


fprintf('   Simulating the ROM in the balanced coordinates...')
T0 = tic; [t, Z2] = ode23s(@(t, z) Fbal(z) + Gbal(z)*u(t), t, z0, opts);
T3 = toc(T0); fprintf("completed in %2.2f seconds. \n", T3)

%% Compute and plot the outputs
y1 = zeros(p,length(t));
y2 = zeros(p,length(t));
for i=1:length(t)
    y1(:,i) = kronPolyEval(h, X1(i,:).');
end
for i=1:length(t)
    y2(:,i) = kronPolyEval(hbal, Z2(i,:).');
end

if r == n
    modelName = 'FOM';
else
    modelName = 'ROM';
end

if nvp.plot
    figure(figNum)
    if degree == 2
        close
        figure(figNum)
        subplot(2,1,1)
        plot(t,y1(n/2-2,:),'DisplayName','FOM output','LineWidth',3.5)
        hold on
        plot(t,y2(n/2-2,:),'--','DisplayName',sprintf('%s output w/ linear transformation',modelName),'LineWidth',3.5)
        
        subplot(2,1,2)
        plot(t,y1(n/2-1,:),'DisplayName','FOM output','LineWidth',3.5)
        hold on
        plot(t,y2(n/2-1,:),'--','DisplayName',sprintf('%s output w/ linear transformation',modelName),'LineWidth',3.5)
    elseif degree == 3
        subplot(2,1,1)
        plot(t,y2(n/2-2,:),':o','MarkerIndices',1:20:length(t),'MarkerSize',6,'DisplayName',sprintf('%s output w/ degree %i transformation',modelName,degree-1))
        
        subplot(2,1,2)
        plot(t,y2(n/2-1,:),':o','MarkerIndices',1:20:length(t),'MarkerSize',6,'DisplayName',sprintf('%s output w/ degree %i transformation',modelName,degree-1))
    else
        subplot(2,1,1)
        plot(t,y2(n/2-2,:),'--+','MarkerIndices',1:15:length(t),'MarkerSize',6,'DisplayName',sprintf('%s output w/ degree %i transformation',modelName,degree-1))
        xlabel('Time, t'); ylabel('$y_1(t)$'); title('Beam tip horizontal displacement')
        legend('Location','southeast'); %ylim(1e-7*[-.35 .05])
        grid on
        
        subplot(2,1,2)
        plot(t,y2(n/2-1,:),'--+','MarkerIndices',1:15:length(t),'MarkerSize',6,'DisplayName',sprintf('%s output w/ degree %i transformation',modelName,degree-1))
        xlabel('Time, t'); ylabel('$y_2(t)$'); title('Beam tip vertical displacement ')
        grid on
        
        set(gcf,"Position",[545 269 774 420])
        exportgraphics(gcf, sprintf('plots/example6_n%i_r%i_d%i_U0%i_y.pdf',n,r,degree,U0*100), 'ContentType', 'vector');
    end
end

error = norm(interp1(t, y2(n/2-2,:), t) - interp1(t, y1(n/2-2,:), t));
fprintf('\n   The output error is: ||yᵣ(t)-y(t)||₂ = %e \n\n', error)

end

