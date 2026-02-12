function [T1, T2, T3, x0] = runExample6_balancedReduction(x0, degree, numEls, r,  U0, verbose, figNum, nvp)
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
    numEls = 3
    r = 4
    U0 = 2e4
    verbose = false
    figNum = 21414836
    nvp.plot = true
    nvp.x0init = zeros(6*numEls,1)
end
clear n nn F2i F2j F2v F3i F3j F3v I21 I22 I31 I32 I33
set(groot,'defaultLineLineWidth',2,'defaultTextInterpreter','latex')
vec = @(x) x(:);

numNodes = numEls + 1; n = 6*numEls;
fprintf('Running Example 6, n=%d, degree %d...\n',n,degree)

%% Get dynamics and define control problem
[f,g,h] = getSystem6_sparse(numEls);
g_artificial = g; h_artificial = h;
g = {g{1}(:,end-1)}; h = {h{1}([n/2-2,n/2-1],:)};

m = size(g{1},2); p = size(h{1},1);
% F = @(x) kronPolyEval(f, x); G = @(x) kronPolyEval(g, x, scenario='G(x)');
F = @(x) sparseKronPolyEval(f, x); G = @(x) g{1};

if isempty(x0)
    % Obtain steady-state Newton iteration for equilibrium point
    fsymmetric = f;
    for i=2:length(fsymmetric)
        fsymmetric{i} = kronMonomialSymmetrize(fsymmetric{i},n,i);
    end
    
    u = [zeros(m-2,1); U0; 0];
    if length(nvp.x0init) ~= n % interpolate lower-order solution as initial guess
        nvp.x0init = vec([zeros(3,2); reshape(nvp.x0init,[],2)]); % add fixed node dofs
        nvp.x0init = vec(interp1(linspace(0,1,length(nvp.x0init)/6), reshape(nvp.x0init,[],6), linspace(0,1,numEls+1), 'linear'));
        nvp.x0init = reshape(nvp.x0init,[],2); % remove fixed node dofs
        nvp.x0init = vec(nvp.x0init(4:end,:));
    end
    x0 = newtonIteration(-g{1}*u, @(x) kronPolyEval(f, x), @(x) sparsejcbn(fsymmetric, x), maxIter=10, z0=nvp.x0init);
end

if verbose
    fprintf('  - The full-order nonlinear model is:\n')
    dispPolyDynamics(f,g,h)
end

%% Compute balanced realization
if verbose
    [fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g_artificial,h_artificial,eta=0,degree=3,transformationDegree=degree-1,f=f,g=g,h=h);
    fprintf('  - The full-order balanced realization for the nonlinear model is:\n')
    dispPolyDynamics(fbal,gbal,hbal,variable='z')
end

tic; fprintf('   Computing balanced reduced-order model...')
[fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g_artificial,h_artificial,r=r,eta=0,degree=3,transformationDegree=degree-1,f=f,g=g,h=h);
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
u = @(t) 0;

%% Simulate original and transformed systems; same input should give same output
% Solve for z0 initial condition
z0 = kronPolyEval(TbalInv,x0);
% z0 = newtonIteration(x0, @(z) kronPolyEval(Tbal,z), @(z) jcbn(Tbal,z),maxIter=100,verbose=true,z0=kronPolyEval(TbalInv,x0));

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
% Simulate both systems
clear nn F2i F2j F2v F3i F3j F3v I21 I22 I31 I32 I33
global T0;
opts = odeset(OutputFcn=@odeprog);
t = [0:.0001:.05];

if n < 90
    fprintf('   Simulating the FOM in the original coordinates...')
    T0 = tic; [t, X1] = ode23s(@(t, x) F(x) + G(x)*u(t), t, x0, opts);
    T2 = toc(T0); fprintf("completed in %2.2f seconds. \n", T2)
else
    T2 = 0;
end

fprintf('   Simulating the ROM in the balanced coordinates...')
T0 = tic; [t, Z2] = ode23s(@(t, z) Fbal(z) + Gbal(z)*u(t), t, z0, opts);
T3 = toc(T0); fprintf("completed in %2.2f seconds. \n", T3)

%% Compute and plot the outputs
if nvp.plot
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
    
    
    figure(figNum); 
    if degree == 2
        close
        figure(figNum); tiledlayout(2,1, 'Padding', 'none', 'TileSpacing', 'compact'); 
        nexttile(1)
        plot(t,y1(1,:),'DisplayName','FOM output','LineWidth',3.5)
        hold on
        plot(t,y2(1,:),'--','DisplayName',sprintf('%s output w/ linear transformation',modelName),'LineWidth',3.5)
        
        nexttile(2)
        plot(t,y1(2,:),'DisplayName','FOM output','LineWidth',3.5)
        hold on
        plot(t,y2(2,:),'--','DisplayName',sprintf('%s output w/ linear transformation',modelName),'LineWidth',3.5)
    elseif degree == 3
        nexttile(1)
        plot(t,y2(1,:),':o','MarkerIndices',1:20:length(t),'MarkerSize',6,'DisplayName',sprintf('%s output w/ degree %i transformation',modelName,degree-1))
        
        nexttile(2)
        plot(t,y2(2,:),':o','MarkerIndices',1:20:length(t),'MarkerSize',6,'DisplayName',sprintf('%s output w/ degree %i transformation',modelName,degree-1))
    else
        nexttile(1)
        plot(t,y2(1,:),'--+','MarkerIndices',1:15:length(t),'MarkerSize',6,'DisplayName',sprintf('%s output w/ degree %i transformation',modelName,degree-1))
        ylabel('$y_1(t)$'); %xlabel('Time, t'); title('Beam tip horizontal displacement')
        legend('Location','southeast'); %ylim(1e-7*[-.35 .05])
        grid on
        
        nexttile(2)
        plot(t,y2(2,:),'--+','MarkerIndices',1:15:length(t),'MarkerSize',6,'DisplayName',sprintf('%s output w/ degree %i transformation',modelName,degree-1))
        ylabel('$y_2(t)$'); xlabel('Time, t'); %title('Beam tip vertical displacement ')
        grid on
        
        set(gcf,"Position",[545 269 774*.75 420*1])
        exportgraphics(gcf, sprintf('plots/example6_n%i_r%i_d%i_U0%i_y.pdf',n,r,degree,U0*100), 'ContentType', 'vector');
    end
end


end

function J = sparsejcbn(F, x)
%jcbn Return the Jacobian J(x) = ∂f(x)/∂x of the function f(x) evaluated at x.
%              f(x) = F₁x + F₂(x⊗x) + ... + Fd(x...⊗x)
%   the Jacobian is given by
%       J(x) = ∂f(x)/∂x = F₁ + 2F₂(I⊗x) + ... + d Fd(I...⊗x)
%%
n = size(x, 1);
% k=1 term
J = full(double(F{1}));
if isempty(find(x,1)) % if x is zero, only the constant term remains, so we are done
    return;
end
% k=2 term, need to iterate over n rows to apply kron-vec identity
for j = 1:n
    if isempty(find(F{2}(j,:),1)); continue; end
    J(j,:) = J(j,:) + 2 * x.' * reshape(F{2}(j,:),n,[]);
end

% k=3 term, need to iterate over n rows to apply kron-vec identity
persistent nn F3i F3j F3v I1 I2
if isempty(nn) || nn ~= length(x)
    nn = length(x);
    for j = 1:nn
        if isempty(find(F{3}(j,:),1)); continue; end
        [F3i{j}, F3j{j}, F3v{j}] = find(reshape(F{3}(j,:),nn,[]));
        [I1{j}, I2{j}] = ind2sub([nn nn], F3j{j});
    end

end
for j = 1:n
    % if isempty(find(F{3}(j,:),1)); continue; end
    if isempty(find(F{3}(j,:),1)); continue; end
    xprod = x(I1{j}) .* x(I2{j});
    J(j,:) = J(j,:) + 3 * accumarray(F3i{j}, F3v{j} .* xprod, [n, 1]).';
end

end


function [x] = sparseKronPolyEval(f,z)
%sparseKronPolyEval Evaluate a Kronecker polynomial with sparse optimization
%   - Assumes f{2} is zero, f{1} and f{3} are only other coefficients
%   - Assumes f{3} is sparse and avoids forming kron(z,z,z)
% Output:   x = f{1}*z + f{3}*(z⊗z⊗z)
%%

% Evaluate linear and quadratic terms normally
x = f{1}*z;

% Use persistent variables to only compute sparsity pattern once
persistent nn F2i F2j F2v F3i F3j F3v I21 I22 I31 I32 I33
if isempty(nn) || nn ~= length(z)
    nn = length(z);
    [F2i, F2j, F2v] = find(f{2});
    [I21, I22] = ind2sub([nn nn], F2j);

    [F3i, F3j, F3v] = find(f{3});
    [I31, I32, I33] = ind2sub([nn nn nn], F3j);
end

% Efficient sparse evaluation of f{2}*(z⊗z)
zprod = z(I21) .* z(I22);
x = x + accumarray(F2i, F2v .* zprod, size(x));

% Efficient sparse evaluation of f{3}*(z⊗z⊗z)
zprod = z(I31) .* z(I32) .* z(I33);
x = x + accumarray(F3i, F3v .* zprod, size(x));

end