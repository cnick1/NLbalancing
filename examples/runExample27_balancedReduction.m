function runExample27_balancedReduction(degree)
%runExample27_balancingTransformation Runs the 3D academic example to
%visualize the nonlinear balanced reduction.
%
%   Usage:  runExample27_balancingTransformation(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description:
%
%
%   References:
%
%   Part of the NLbalancing repository.
%%
% close all;
arguments
    degree = 4
end
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 27\n')

%% Get system dynamics
% Get dynamics
n = 12; r=3;
[f, g, h, xg] = getSystem27(n-1,.75,1,-1);
% f = f(1); g = g(1); h = h(1);

% Get value function/controller
q = 0.5; R = 1;
[~, K] = ppr(f, g, q, R, 2);
f{1} = f{1}+g{1}*K{1};

%% Compute balanced realization
[fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g,h,r=r,transformationDegree=degree-1);
TbalInv = transformationInverse(Tbal);

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
Fbal = @(z) kronPolyEval(fbal, z);

%% Simulate original and transformed systems; same input should give same output
% Solve for z0 initial condition
x0 = 1 - 12*xg.^2 + 12* xg.^3 - 1*xg.^12; x0 = 1*x0.';

z0 = kronPolyEval(TbalInv,x0,1); % Actually, the Newton iteration approach could never work :)
fprintf('         -> Initial condition error: %2.2e \n', norm(kronPolyEval(Tbal,z0,1)-x0))

%% Apply reduction by eliminating z3 (set it and its derivative to zero)

% Simulate both systems
global T0; opts = odeset(OutputFcn=@odeprog); 
T0 = tic; [t1, X1] = ode45(@(t, x) F(x), [0, 5], x0,opts);
T0 = tic; [t2, Z2] = ode45(@(t, z) Fbal(z), [0, 5], z0,opts);

% Convert z trajectories to x coordinates
X2 = zeros(length(t2),n);
for i=1:length(t2)
    X2(i,:) = kronPolyEval(Tbal,Z2(i,:).');
end

% Convert x trajectories to z coordinates
Z1 = zeros(length(t1),r);
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

% Plot outputs
hold on;
plot(t1,y1)
plot(t2,y2,'--')
plot(t2,y2,':')

title('Model output'); xlabel('Time t')
ylabel('y(t)')
legend('FOM output','ROM output using hbar(z)','ROM output h(x)')

fprintf('The output error is: %f \n', norm(interp1(t2, y2, 0:.1:5) - interp1(t1, y1, 0:.1:5)))


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