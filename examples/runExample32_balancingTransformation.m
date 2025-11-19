function runExample32_balancingTransformation(degree,lim,reduction)
%runExample32_balancingTransformation Runs the 3D academic example to
%visualize the nonlinear balanced reduction.
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
%   controllability (since it drives the other states). We apply the
%   following nonlinear transformation to the LTI balanced realization
%           x = 𝟁(z) = [z₁;   z₂;   z₃ + z₁² + z₂² + z₁³]
%   to obtain nonlinear dynamics that we can test our method on.
%   We compute the energy functions, the input-normal/output-diagonal
%   transformation, and then the true balancing transformation, given by the
%   composition x = ̅Φ(z̄(z̄) = Φ(𝝋(z̄)). Since the transformation was
%   applied to the balanced realization, the balancing transformation will
%   be the inverse of this transformation, which we know analytically to be
%
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

fprintf('Running Example 32, polynomial balanced realization...\n')

% x0 = [1 -1 -3].'*lim;
% rng(1);
% x0 = [rand(1);rand(1);0];
x0 = [-1;-2;0];
x0 = [x0(1);x0(2);-x0(1)^2-x0(1)^3-x0(2)^2];

%% Get system dynamics
[f, g, h] = getSystem32(transform=true);

%% Compute balanced realization
[fbal,gbal,hbal,Tbal,sigmaSquared] = getBalancedRealization(f,g,h,eta=0,degree=3,transformationDegree=degree-1);
TbalInv = transformationInverse(Tbal);

fprintf('  - The balanced realization for the nonlinear model is:\n')
dispPolyDynamics(fbal,gbal,hbal)

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
Fbal = @(z) kronPolyEval(fbal, z);

%% Simulate original and transformed systems; same input should give same output
% Solve for z0 initial condition
z0 = kronPolyEval(TbalInv,x0);
fprintf(['   Simulating the system in the original vs transformed coordinates for initial condition x0 = [', repmat('%2.2e ', 1, numel(x0)), '], ...\n'], x0)
fprintf(['   ... the transformed initial condition is z0 = [', repmat('%2.2e ', 1, numel(z0)), '] '], z0)
fprintf(' (error: %2.2e). \n', norm(kronPolyEval(Tbal,z0)-x0))

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

% Convert x trajectories to z coordinates
Z1 = zeros(size(X1));
for i=1:length(t1)
    Z1(i,:) = kronPolyEval(TbalInv,X1(i,:).');
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
y3 = zeros(size(t2));
for i=1:length(t2)
    y3(i) = kronPolyEval(h,X2(i,:).');
end

% Plot state trajectories
figure('Position', [827 220 560 420]);
subplot(3,1,1); hold on;
plot(t1,X1)
plot(t2,X2,'--')

title('Reconstructed state trajectories'); xlabel('Time t')
ylabel('x_i(t)')
legend('FOM x_1','FOM x_2','FOM x_3','ROM x_1','ROM x_2','ROM x_3')

% Plot transformed state trajectories
subplot(3,1,2); hold on;
plot(t1,Z1)
plot(t2,Z2,'--')

title('Reconstructed state trajectories'); xlabel('Time t')
ylabel('z_i(t)')
legend('FOM z_1','FOM z_2','FOM z_3','ROM z_1','ROM z_2','ROM z_3')

% Plot outputs
subplot(3,1,3); hold on;
plot(t1,y1)
plot(t2,y2,'--')
plot(t2,y2,':')

title('Model output'); xlabel('Time t')
ylabel('y(t)')
legend('FOM output','ROM output using hbar(z)','ROM output h(x)')

fprintf('\n   The output error is: ||yᵣ(t)-y(t)||₂ = %f \n\n', norm(interp1(t2, y2, 0:.1:5) - interp1(t1, y1, 0:.1:5)))

%% Manifold figure plot
dx = 0.05; lim = 2;
[z1,z2,z3] = meshgrid(-lim:dx:lim,-lim:dx:lim,0);
x1 = zeros(size(z1)); x2 = x1; x3 = x1;
for i=1:numel(z1)
    [x1(i), x2(i), x3(i)] = kronPolyEval(Tbal,[z1(i);z2(i);z3(i)]);
end

% figure;
% surf(z1,z2,z3)
figure;
surf(x1,x2,x3)
xlim([-2 3])
ylim([-2 2])
zlim([-20 5])
drawnow
% return

    hold on;
    plot3(X1(:,1),X1(:,2),X1(:,3),'g',DisplayName='FOM solution')
    plot3(X2(:,1),X2(:,2),X2(:,3),'r',DisplayName='ROM solution')

if degree == 2
    fprintf('    -> The figure confirms that the reduced-order solution trajectories on the linear balanced subspace fail to capture the full dynamics. \n\n')
else
    fprintf('    -> The figure confirms that the reduced-order solution trajectories on the balanced manifold produce the same output as the original dynamics. \n\n')
end
%% Simulate the system's response to a random noise input
% Instead of simulating the nonlinear system, it is much faster to simulate
% the linear system and then transform the solution
if false && reduction
    [fl, gl, hl] = getSystem32(transform=false);
    [fl, gl, hl] = getBalancedRealization(fl,gl,hl,eta=0,degree=1);
    
    global T0; opts = odeset(OutputFcn=@odeprog); T0 = tic;
    [T, Z] = ode45(@(t, z) (fl{1}*z + gl{1}*randn(1)), [0, 10], [0;0;0],opts);
    X = zeros(size(Z)); Z = 100*Z; % Scale the output to better see the manifold
    for i=1:length(T)
        X(i,:) = [Z(i,1), Z(i,2), Z(i,3) - Z(i,1)^2 - Z(i,2)^2 - Z(i,1)^3]; % Apply the inverse transformation to get the nonlinear solution
    end
    
    hold on;
    plot3(X(:,1),X(:,2),X(:,3),'r')
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3')
    title('Balanced manifold and white noise response')
end
end

function status = odeprog(t, y, flag)
% ODEPROG Custom progress bar for ode solver
% Use with odeset: opts = odeset('OutputFcn',@odeprog);
persistent T1 nSteps lastPct lastUpdateTime
global T0

status = false;

switch flag
    case 'init'
        % Initialize progress bar
        elapsed = toc(T0);
        T1 = t(end);
        nSteps = 50; % Number of blocks in the progress bar
        lastPct = -1;
        lastUpdateTime = 0;
        fprintf(' |%s|  (elapsed: %5i s, remaining: ----- s)', repmat(' ',1,nSteps), round(elapsed));
        
    case ''
        % ODE solver step
        if isempty(t), return; end
        tNow = t(end);
        pct = min(100, max(0, 100 * tNow / T1));
        block = floor(pct / (100/nSteps));
        elapsed = toc(T0);
        eta = (elapsed / max(tNow,eps)) * (T1 - tNow); % avoid divide-by-zero
        needsUpdate = pct - lastPct >= 2 || block == nSteps;
        timeSinceLast = elapsed - lastUpdateTime;
        
        if needsUpdate || timeSinceLast >= 1
            bar = [repmat('-',1,block), repmat(' ',1,nSteps-block)];
            fprintf(repmat('\b',1,93));
            fprintf(' |%s|  (elapsed: %5i s, remaining: %5i s)', bar, round(elapsed), min(round(eta),99999));
            if needsUpdate
                lastPct = pct;
            end
            lastUpdateTime = elapsed;
        end
        
    case 'done'
        % Finalize
        % elapsed = toc(T0);
        % bar = repmat('-',1,nSteps);
        fprintf(repmat('\b',1,93));
        % fprintf(' |%s|  (elapsed: %5i s, remaining:     0 s)\n', bar, round(elapsed));
        clear T0 T1 nSteps lastPct lastUpdateTime
end
end