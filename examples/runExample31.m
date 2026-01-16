function runExample31(degree,lim)
%runExample31 Runs 2D pendulum example to compare with Newman's results
%
%   Usage:  runExample31(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description: We consider the pendulum as described in [1-3] by
%
%             ẍ + G/L sin(x) + k / m L² x + b/mL² ẋ = 1/mL² u(t)
%
%   This can be put in first-order form with x₁ = x, x₂ = ẋ as
%
%             ẋ₁ = x₂
%      (1)    ẋ₂ = -G/L sin(x₁) - k / m L² x₁ - b/mL² x₂ + u(t)
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
%   References: [1] A. J. Newman and P. S. Krishnaprasad, “Computation for
%                   nonlinear balancing,” University of Maryland, College
%                   Park, 1998.
%               [2] A. J. Newman and P. S. Krishnaprasad, “Computing
%                   balanced realizations for nonlinear systems,” University
%                   of Maryland, College Park, 2000.
%               [3] A. J. Newman, “Modeling and reduction with applications
%                   to semiconductor processing,” University of Maryland,
%                   College Park, 1999.
%
%   Part of the NLbalancing repository.
%%
clc; close all;
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex','defaultAxesTickLabelInterpreter', 'latex', 'defaulttextinterpreter', 'latex', 'defaultLegendInterpreter', 'latex');

fprintf('Running Example 31\n')

if nargin < 2
    lim = 1;
    if nargin < 1
        degree = 10;
    end
end

%% Get system dynamics
[f, g, h, FofXU] = getSystem31(degree-1);

%%  Compute the energy functions
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~ Computing energy functions:  ~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
eta = 0;
[v] = approxPastEnergy(f, g, h, eta=eta, degree=degree);
[w, K] = approxFutureEnergy(f, g, h, eta=eta, degree=degree);

% Plot the energy functions
N = 31;
xPlot = linspace(-1, 1, N); yPlot = linspace(-1, 1, N);
[X, Y] = meshgrid(xPlot, yPlot);

ePast = zeros(N, N); eFuture = zeros(N, N);
for i = 1:N
    for j = 1:N
        x = [X(i, j); Y(i, j)];
        ePast(i, j) = 0.5 * kronPolyEval(v, x, degree=degree);
        eFuture(i, j) = 0.5 * kronPolyEval(w, x, degree=degree);
    end
end
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex','defaultAxesTickLabelInterpreter', 'latex', 'defaulttextinterpreter', 'latex', 'defaultLegendInterpreter', 'latex');

figure('Position',[4.3333 57 362.6667 360]);
% contourf(X, Y, ePast, 16, 'w'); %colormap(flip(YlGnBuRescaled))
surfc(X, Y, ePast); %colormap(flip(YlGnBuRescaled))
xlabel('$x_1$'); ylabel('$x_2$');  title('Controllability Energy')
% xticks(-1:1); yticks(-1:1)
set(gca, 'FontSize', 16)
% h = colorbar('FontSize', 16, 'Location','southoutside');
% clim([0 80]); set(h, 'ylim', [0 80])

figure('Position',[365.6667 57 362.6667 360]);
% contourf(X, Y, eFuture, 16, 'w'); %colormap(flip(YlGnBuRescaled))
surfc(X, Y, eFuture); %colormap(flip(YlGnBuRescaled))
xlabel('$x_1$'); ylabel('$x_2$');  title('Observability Energy')
% xticks(-1:1); yticks(-1:1)
set(gca, 'FontSize', 16)
% h = colorbar('FontSize', 16,'Location','southoutside');
% clim([0 1.5]); set(h, 'ylim', [0 1.5])

% Plot the past and future HJB residuals
xPlot = linspace(-1, 1, N); yPlot = linspace(-1, 1, N);
[X, Y] = meshgrid(xPlot, yPlot);

vRES = computeResidualPastHJB(f, g, h, eta, v, degree, 1, N);
wRES = computeResidualFutureHJB(f, g, h, eta, w, degree, 1, N);

figure('Position',[880 287 363.3333 360]);
% contourf(X, Y, abs(vRES), 16, 'w'); %colormap(flip(YlGnBuRescaled))
surf(X, Y, vRES);
xlabel('$x_1$'); ylabel('$x_2$'); title('Controllability Energy Residual')
zlim([-.8 .6])
% xticks(-1:1); yticks(-1:1)
% h = colorbar('FontSize', 16,'Location','southoutside');
set(gca, 'FontSize', 16)
% clim([0,2.5e-11])
% set(h,'YTick',0:1e-11:2e-11,'TickLabels',{'0','$1\cdot10^{-11}$','$2\cdot10^{-11}$'})
fprintf('\nThe residual of the HJB equation on the unit square is %g\n', norm(vRES, 'inf'));

figure('Position',[1270 287 363.3333 360]);
% pcolor(X, Y, abs(wRES)); shading interp; %colormap(flip(YlGnBuRescaled))
surf(X, Y, wRES);
xlabel('$x_1$'); ylabel('$x_2$'); title('Observability Energy Residual')
zlim([-.3 .4])
% xticks(-1:1); yticks(-1:1)
% h = colorbar('FontSize', 16,'Location','southoutside');
set(gca, 'FontSize', 16)
% set(h,'YTick',0:2e-16:4e-16,'TickLabels',{'0','$2\cdot10^{-16}$','$4\cdot10^{-16}$'})
fprintf('\nThe residual of the HJB equation on the unit square is %g\n', norm(wRES, 'inf'));
drawnow



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
figure("Position", [731 57 997.3333 420]);
subplot(1,2,2); title(sprintf("z for z ∈[-%1.1f,%1.1f] $\\times$ [-%1.1f,%1.1f] ",lim,lim,lim,lim),'Interpreter','none')
hold on; axis equal; xlabel("$z_1$"); ylabel("$z_2$")
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i))
    plot(xV(:,i),yV(:,i))
end
% xlim([-.5 .5]);ylim([-.5 .5])

subplot(1,2,1); title(sprintf("x = Φ(z) for z ∈[-%1.1f,%1.1f] $\\times$ [-%1.1f,%1.1f] ",lim,lim,lim,lim),'Interpreter','none')
hold on; %axis equal;
xlabel("$x_1$"); ylabel("$x_2$")
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i))
    plot(xVtr(:,i),yVtr(:,i))
end
% xlim([-.2 .45]);ylim([-.2 .45])

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
uPPR = @(x) kronPolyEval(K, x);
% FofXU = @(x, u) [x(2); -sin(x(1))] + g{1} * u;
% FofXU = @(x) [x(2); -sin(x(1))];
FtofZU = @(z,u) PhiBarJacobian(z,TinOd,sigmaSquared)\FofXU(PhiBar(z,TinOd,sigmaSquared),u);

% theta0 = -pi;
% x0 = [theta0; 0];
x0 = [0.1; 0.1];

% Solve for z0 initial condition with a Newton type iteration
z0 = newtonIteration(x0, @(z) PhiBar(z,TinOd,sigmaSquared), @(z) PhiBarJacobian(z,TinOd,sigmaSquared),maxIter=10,verbose=true);

%% Simulate both systems, u = 0
[t, X1] = ode45(@(t, x) FofXU(x,uPPR(x)), [0, 40], x0);
[~, Z] = ode45(@(t, z) FtofZU(z,uPPR(x)), t, z0);

subplot(1,2,1)
plot(X1(:,1),X1(:,2),'g','LineWidth',1.5)
subplot(1,2,2)
plot(Z(:,1),Z(:,2),'r--','LineWidth',1.5)


% Transform Z trajectory into X coordinates to compare
X2 = zeros(size(Z));
for i = 1:length(Z)
    X2(i,:) = PhiBar(Z(i,:).',TinOd,sigmaSquared);
end

subplot(1,2,1)
plot(X2(:,1),X2(:,2),'r--','LineWidth',1.5)
xlim([-pi pi])
ylim([-pi pi])
drawnow

figure
plot(t,X1(:,1))
hold on
plot(t,X2(:,1),'--')
grid on
xlabel("Time"); ylabel("Position")

%% Simulate both systems, u = 0.5 sin(t/π)
[t, X1] = ode45(@(t, x) FofXU(x, 0.5*sin(t/pi)), [0, 40], x0);
[~, Z] = ode45(@(t, z) FtofZU(z, 0.5*sin(t/pi)), t, z0);

figure("Position", [731 57 997.3333 420]);
subplot(1,2,2); title(sprintf("z for z ∈[-%1.1f,%1.1f] $\\times$ [-%1.1f,%1.1f] ",lim,lim,lim,lim),'Interpreter','none')
hold on; axis equal; xlabel("$z_1$"); ylabel("$z_2$")
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i))
    plot(xV(:,i),yV(:,i))
end
% xlim([-.5 .5]);ylim([-.5 .5])

subplot(1,2,1); title(sprintf("x = Φ(z) for z ∈[-%1.1f,%1.1f] $\\times$ [-%1.1f,%1.1f] ",lim,lim,lim,lim),'Interpreter','none')
hold on; %axis equal;
xlabel("$x_1$"); ylabel("$x_2$")
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i))
    plot(xVtr(:,i),yVtr(:,i))
end

subplot(1,2,1)
plot(X1(:,1),X1(:,2),'g','LineWidth',1.5)
subplot(1,2,2)
plot(Z(:,1),Z(:,2),'r--','LineWidth',1.5)


% Transform Z trajectory into X coordinates to compare
X2 = zeros(size(Z));
for i = 1:length(Z)
    X2(i,:) = PhiBar(Z(i,:).',TinOd,sigmaSquared);
end

subplot(1,2,1)
plot(X2(:,1),X2(:,2),'r--','LineWidth',1.5)
% xlim([-pi pi]); ylim([-pi pi])
xlim([-.1 .1]); ylim([-.1 .1])
drawnow

figure
plot(t,X1(:,1))
hold on
plot(t,X2(:,1),'--')
grid on
xlabel("Time"); ylabel("Position")


end

