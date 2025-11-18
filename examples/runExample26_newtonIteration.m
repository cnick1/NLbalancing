function runExample26_newtonIteration(degree,lim)
%runExample26_newtonIteration Runs 2D inverted pendulum example to
%   visualize the nonlinear balancing transformations.
%
%   Usage:  runExample26_newtonIteration(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description: The basic equation for the inverted pendulum is
%
%               ẍ - sin(x) = u(t)
%
%   This can be put in first-order form with x₁ = x, x₂ = ẋ as
%
%             ẋ₁ = x₂
%      (1)    ẋ₂ = sin(x₁) + u(t)
%              y = x₁
%
%   and sin(x₁) can be approximated as
%
%             sin(x₁) = x₁ - x₁³/6 + x₁⁵/120 - x₁⁷/5040 + x₁⁹/362880 + ...
%
%   We compute the energy functions, the input-normal/output-diagonal
%   transformation, and then the true balancing transformation, given by the
%   composition x = ̅Φ(z̄(z̄) = Φ(𝝋(z̄)). The mapping from the z̄ coordinates
%   to the x coordinates is visualized by forming a grid in the z̄ coordinates
%   and mapping that grid to the x coordinates.
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
% close all;
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 26\n')

if nargin < 2
    lim = 5;
    if nargin < 1
        degree = 4;
    end
end

%% Get system dynamics
[f, g, h] = getSystem26(degree-1);
for i=1:length(f) % make it the stable pendulum
    f{i}(2,:) = -f{i}(2,:);
end

%%  Compute the energy functions
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~ Computing energy functions:  ~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
eta = 0.00001;
[v] = approxPastEnergy(f, g, h, eta=eta, degree=degree);
[w, K] = approxFutureEnergy(f, g, h, eta=eta, degree=degree);

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
fprintf("\n ~~~~~~~~~~~ Computing transformation and singular value functions:  ~~~~~~~~~~~~ \n")
[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree=degree-1, verbose=true);

%% Plot grid transformations
% Parameters
numLines = 41; numPoints = 201;

% Generate original z coordinates
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

% Compute transformed coordinates
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
for i=1:length(xH(:))
    % temp = kronPolyEval(TinOd, [xH(i);yH(i)]);
    temp = PhiBar([xH(i);yH(i)],TinOd,sigmaSquared);
    xHtr(i) = temp(1); yHtr(i) = temp(2);
    
    % temp = kronPolyEval(TinOd, [xV(i);yV(i)]);
    temp = PhiBar([xV(i);yV(i)],TinOd,sigmaSquared);
    xVtr(i) = temp(1); yVtr(i) = temp(2);
end

% Prepare figure
figure("Position", [185 337.6667 997.3333 420]);
subplot(1,2,2); title(sprintf("z for z ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
hold on; axis equal; xlabel("z_1"); ylabel("z_2")
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i))
    plot(xV(:,i),yV(:,i))
end
% xlim([-.5 .5]);ylim([-.5 .5])

subplot(1,2,1); title(sprintf("x = Φ(z) for z ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
hold on; %axis equal;
xlabel("x_1"); ylabel("x_2")
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i))
    plot(xVtr(:,i),yVtr(:,i))
end
% xlim([-.2 .45]);ylim([-.2 .45])

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
uPPR = @(x) kronPolyEval(K, x);
FofXU = @(x) [x(2); -sin(x(1))] + g{1} * uPPR(x);
% FofXU = @(x) [x(2); -sin(x(1))];
FtofZU = @(z) PhiBarJacobian(z,TinOd,sigmaSquared)\FofXU(PhiBar(z,TinOd,sigmaSquared));

theta0 = -pi;
x0 = [theta0; 0];

% Solve for z0 initial condition with a Newton type iteration
z0 = newtonIteration(x0, @(z) PhiBar(z,TinOd,sigmaSquared), @(z) PhiBarJacobian(z,TinOd,sigmaSquared),maxIter=10,verbose=true);

% Simulate both systems
[~, X1] = ode45(@(t, x) FofXU(x), [0, 50], x0);
[~, Z] = ode45(@(t, z) FtofZU(z), [0, 50], z0);

subplot(1,2,1)
plot(X1(:,1),X1(:,2),'g','LineWidth',1.5)
subplot(1,2,2)
plot(Z(:,1),Z(:,2),'r--','LineWidth',1.5)


%% Transform Z trajectory into X coordinates to compare
X2 = zeros(size(Z));
for i = 1:length(Z)
    X2(i,:) = PhiBar(Z(i,:).',TinOd,sigmaSquared);
end

subplot(1,2,1)
plot(X2(:,1),X2(:,2),'r--','LineWidth',1.5)
xlim([-pi pi])
ylim([-pi pi])
drawnow
end


