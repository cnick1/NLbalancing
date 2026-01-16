function runExample31_nonlinearBalancing()
%runExample31_nonlinearBalancing Runs 2D pendulum example to compare with Newman's results
%
%   Usage:  runExample31_nonlinearBalancing()
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

fprintf('Running Example 31\n')

%% Get system dynamics
degree = 4;
[f, g, h, FofXU] = getSystem31(degree-1);

%%  Compute the energy functions
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

figure('Position',[4.3333 57 387 318]);
surfc(X, Y, ePast); 
xlabel('$x_1$'); ylabel('$x_2$'); title('Controllability Energy')

figure('Position',[365.6667 57 387 318]);
surfc(X, Y, eFuture);
xlabel('$x_1$'); ylabel('$x_2$'); title('Observability Energy')

% Plot the past and future HJB residuals
xPlot = linspace(-1, 1, N); yPlot = linspace(-1, 1, N);
[X, Y] = meshgrid(xPlot, yPlot);

vRES = computeResidualPastHJB(f, g, h, eta, v, degree, 1, N);
wRES = computeResidualFutureHJB(f, g, h, eta, w, degree, 1, N);

figure('Position',[880 287 387 318]);
contourf(X, Y, vRES); colorbar; axis equal
xlabel('$x_1$'); ylabel('$x_2$'); title('Controllability Energy Residual')
zlim([-.8 .6])

figure('Position',[1270 287 387 318]);
contourf(X, Y, wRES); colorbar; axis equal
xlabel('$x_1$'); ylabel('$x_2$'); title('Observability Energy Residual')
zlim([-.3 .4])
drawnow

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
[fbal,gbal,hbal,Tbal] = getBalanceThenReduceRealization(f,g,h,eta=0,transformationDegree=degree-1,verbose=true);
TbalInv = transformationInverse(Tbal);
% FofXU = @(x,u) kronPolyEval(f,x) + kronPolyEval(g,x,scenario='G(x)')*u;
FtofZU = @(z,u) kronPolyEval(fbal,z) + kronPolyEval(gbal,z,scenario='G(x)')*u;

%% Simulate dynamics: impulse-type response
% Solve for z0 initial condition with a Newton type iteration
x0 = [0.1; 0.1];
z0 = kronPolyEval(TbalInv,x0);

[t, X1] = ode45(@(t, x) FofXU(x,0), [0, 40], x0);
[~, Z1] = ode45(@(t, z) FtofZU(z,0), t, z0);

% Transform Z trajectory into X coordinates to compare
X2 = zeros(size(Z1));
for i = 1:length(Z1)
    X2(i,:) = kronPolyEval(Tbal,Z1(i,:).');
end
Z2 = zeros(size(X1));
for i = 1:length(X1)
    Z2(i,:) = kronPolyEval(TbalInv,X1(i,:).');
end

% Plot grid transformations
% Parameters
numLines = 27; numPoints = 201;
lim = 1;
% Generate original z coordinates
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

% Compute transformed coordinates
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
zxHtr = zeros(size(xH)); zyHtr = zeros(size(yH));
zxVtr = zeros(size(xV)); zyVtr = zeros(size(yV));
for i=1:length(xH(:))
    [xHtr(i), yHtr(i)] = kronPolyEval(Tbal,[xH(i);yH(i)]);
    [xVtr(i), yVtr(i)] = kronPolyEval(Tbal,[xV(i);yV(i)]);

    [zxHtr(i), zyHtr(i)] = kronPolyEval(TbalInv,[xH(i);yH(i)]);
    [zxVtr(i), zyVtr(i)] = kronPolyEval(TbalInv,[xV(i);yV(i)]);
end

% Prepare figure
% Generate figure
figure('Position',[850 442 396 300]); hold on
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i),'Color', [0.75 0.75 0.75])
    plot(xVtr(:,i),yVtr(:,i),'Color', [0.75 0.75 0.75])
end
set(gca, 'ColorOrderIndex', 1);
plot(X1(:,1),X1(:,2),'LineWidth',1.5)
plot(X2(:,1),X2(:,2),'--','LineWidth',1.5)
axis equal;
xlim([-.2 .2]); ylim([-.2 .2])
xlabel("$x_1$"); ylabel("$x_2$")
exportgraphics(gca, 'plots/example31_impulse_x.pdf', 'ContentType', 'vector');

figure('Position',[1250 442 396 300]); hold on
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i),'Color', [0.75 0.75 0.75])
    plot(xV(:,i),yV(:,i),'Color', [0.75 0.75 0.75])
end
set(gca, 'ColorOrderIndex', 1);
plot(Z2(:,1),Z2(:,2),'LineWidth',1.5)
plot(Z1(:,1),Z1(:,2),'--','LineWidth',1.5)
axis equal;
xlim([-.5 .5]); ylim([-.5 .5])
xlabel("$z_1$"); ylabel("$z_2$")
exportgraphics(gca, 'plots/example31_impulse_z.pdf', 'ContentType', 'vector');

figure('Position',[290 442 560 300]); hold on
plot(t,X1(:,1),'LineWidth',1.5,'DisplayName','Original realization solution')
plot(t,X2(:,1),'--','LineWidth',1.5,'DisplayName','Balanced realization solution')
legend; grid on
xlabel("Time, t"); ylabel("Angular Position, $x_1$")
exportgraphics(gca, 'plots/example31_impulse.pdf', 'ContentType', 'vector');

%% Simulate dynamics: u = 0.5 sin(t/π)
% Solve for z0 initial condition with a Newton type iteration
x0 = [0.1; 0.1];
z0 = kronPolyEval(TbalInv,x0);

[t, X1] = ode45(@(t, x) FofXU(x, 0.5*sin(t/pi)), [0, 40], x0);
[~, Z1] = ode45(@(t, z) FtofZU(z, 0.5*sin(t/pi)), t, z0);

% Transform Z trajectory into X coordinates to compare
X2 = zeros(size(Z1));
for i = 1:length(Z1)
    X2(i,:) = kronPolyEval(Tbal,Z1(i,:).');
end
Z2 = zeros(size(X1));
for i = 1:length(X1)
    Z2(i,:) = kronPolyEval(TbalInv,X1(i,:).');
end

% Plot grid transformations
% Parameters
numLines = 27; numPoints = 201;

% Generate original z coordinates
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

% Compute transformed coordinates
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
zxHtr = zeros(size(xH)); zyHtr = zeros(size(yH));
zxVtr = zeros(size(xV)); zyVtr = zeros(size(yV));
for i=1:length(xH(:))
    [xHtr(i), yHtr(i)] = kronPolyEval(Tbal,[xH(i);yH(i)]);
    [xVtr(i), yVtr(i)] = kronPolyEval(Tbal,[xV(i);yV(i)]);

    [zxHtr(i), zyHtr(i)] = kronPolyEval(TbalInv,[xH(i);yH(i)]);
    [zxVtr(i), zyVtr(i)] = kronPolyEval(TbalInv,[xV(i);yV(i)]);
end

% Prepare figure
% Generate figure
figure('Position',[850 442 396 300]); hold on
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i),'Color', [0.75 0.75 0.75])
    plot(xVtr(:,i),yVtr(:,i),'Color', [0.75 0.75 0.75])
end
set(gca, 'ColorOrderIndex', 1);
plot(X1(:,1),X1(:,2),'LineWidth',1.5)
plot(X2(:,1),X2(:,2),'--','LineWidth',1.5)
axis equal;
xlim([-.2 .2]); ylim([-.2 .2])
xlabel("$x_1$"); ylabel("$x_2$")
exportgraphics(gca, 'plots/example31_sin1_x.pdf', 'ContentType', 'vector');

figure('Position',[1250 442 396 300]); hold on
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i),'Color', [0.75 0.75 0.75])
    plot(xV(:,i),yV(:,i),'Color', [0.75 0.75 0.75])
end
set(gca, 'ColorOrderIndex', 1);
plot(Z2(:,1),Z2(:,2),'LineWidth',1.5)
plot(Z1(:,1),Z1(:,2),'--','LineWidth',1.5)
axis equal;
xlim([-.5 .5]); ylim([-.5 .5])
xlabel("$z_1$"); ylabel("$z_2$")
exportgraphics(gca, 'plots/example31_sin1_z.pdf', 'ContentType', 'vector');

figure('Position',[290 442 560 300]); hold on
plot(t,X1(:,1),'LineWidth',1.5,'DisplayName','Original realization solution')
plot(t,X2(:,1),'--','LineWidth',1.5,'DisplayName','Balanced realization solution')
legend; grid on
xlabel("Time, t"); ylabel("Angular Position, $x_1$")
exportgraphics(gca, 'plots/example31_sin1.pdf', 'ContentType', 'vector');


%% Simulate dynamics: u = 5 sin(t/π)
% Solve for z0 initial condition with a Newton type iteration
x0 = [1; 1];
z0 = kronPolyEval(TbalInv,x0);

[t, X1] = ode45(@(t, x) FofXU(x, 1*sin(t/pi)), [0, 40], x0);
[~, Z1] = ode45(@(t, z) FtofZU(z, 1*sin(t/pi)), t, z0);

% Transform Z trajectory into X coordinates to compare
X2 = zeros(size(Z1));
for i = 1:length(Z1)
    X2(i,:) = kronPolyEval(Tbal,Z1(i,:).');
end
Z2 = zeros(size(X1));
for i = 1:length(X1)
    Z2(i,:) = kronPolyEval(TbalInv,X1(i,:).');
end

% Plot grid transformations
% Parameters
numLines = 27; numPoints = 201;
lim = 6.75;
% Generate original z coordinates
[yH, xH] = meshgrid(linspace(-lim, lim, numLines), linspace(-2*lim, 2*lim, 2*numPoints)); % Vertical lines
[xV, yV] = meshgrid(linspace(-2*lim, 2*lim, 2*numLines), linspace(-lim, lim, numPoints)); % Horizontal lines

% Compute transformed coordinates
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
zxHtr = zeros(size(xH)); zyHtr = zeros(size(yH));
zxVtr = zeros(size(xV)); zyVtr = zeros(size(yV));
for i=1:length(xH(:))
    [xHtr(i), yHtr(i)] = kronPolyEval(Tbal,[xH(i);yH(i)]);
    [xVtr(i), yVtr(i)] = kronPolyEval(Tbal,[xV(i);yV(i)]);
    
    [zxHtr(i), zyHtr(i)] = kronPolyEval(TbalInv,[xH(i);yH(i)]);
    [zxVtr(i), zyVtr(i)] = kronPolyEval(TbalInv,[xV(i);yV(i)]);
end

% Prepare figure
% Generate figure
figure('Position',[850 442 396 300]); hold on
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i),'Color', [0.75 0.75 0.75])
end
for i=1:size(xV,2)
    plot(xVtr(:,i),yVtr(:,i),'Color', [0.75 0.75 0.75])
end
set(gca, 'ColorOrderIndex', 1);
plot(X1(:,1),X1(:,2),'LineWidth',1.5)
plot(X2(:,1),X2(:,2),'--','LineWidth',1.5)
axis equal;
xlim([-1.75 1.75]); ylim([-1.75 1.75])
xlabel("$x_1$"); ylabel("$x_2$")
exportgraphics(gca, sprintf('plots/example31_sin2_x_d%i.pdf',degree), 'ContentType', 'vector');

figure('Position',[1250 442 396 300]); hold on
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i),'Color', [0.75 0.75 0.75])
end
for i=1:size(xV,2)
    plot(xV(:,i),yV(:,i),'Color', [0.75 0.75 0.75])
end
set(gca, 'ColorOrderIndex', 1);
plot(Z2(:,1),Z2(:,2),'LineWidth',1.5)
plot(Z1(:,1),Z1(:,2),'--','LineWidth',1.5)
axis equal;
xlim([-5 5]); ylim([-5 5])
xlabel("$z_1$"); ylabel("$z_2$")
exportgraphics(gca, sprintf('plots/example31_sin2_z_d%i.pdf',degree), 'ContentType', 'vector');

figure('Position',[290 442 560 300]); hold on
plot(t,X1(:,1),'LineWidth',1.5,'DisplayName','Original realization solution')
plot(t,X2(:,1),'--','LineWidth',1.5,'DisplayName','Balanced realization solution')
legend; grid on
xlabel("Time, t"); ylabel("Angular Position, $x_1$")
exportgraphics(gca, sprintf('plots/example31_sin2_d%i.pdf',degree), 'ContentType', 'vector');


% Plot grid transformations
% Parameters
numLines = 27; numPoints = 201; lim = 10;
% Generate original z coordinates
[yH, xH] = meshgrid(linspace(-lim, lim, numLines), linspace(-2*lim, 2*lim, 2*numPoints)); % Vertical lines
[xV, yV] = meshgrid(linspace(-2*lim, 2*lim, 2*numLines), linspace(-lim, lim, numPoints)); % Horizontal lines

% Compute transformed coordinates
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
zxHtr = zeros(size(xH)); zyHtr = zeros(size(yH));
zxVtr = zeros(size(xV)); zyVtr = zeros(size(yV));
for i=1:length(xH(:))
    [xHtr(i), yHtr(i)] = kronPolyEval(Tbal,[xH(i);yH(i)]);
    [xVtr(i), yVtr(i)] = kronPolyEval(Tbal,[xV(i);yV(i)]);
    
    [zxHtr(i), zyHtr(i)] = kronPolyEval(TbalInv,[xH(i);yH(i)]);
    [zxVtr(i), zyVtr(i)] = kronPolyEval(TbalInv,[xV(i);yV(i)]);
end

% Prepare figure
% Generate figure
figure('Position',[850 442 396 300]); hold on
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i),'Color', [0.75 0.75 0.75])
end
for i=1:size(xV,2)
    plot(xVtr(:,i),yVtr(:,i),'Color', [0.75 0.75 0.75])
end
set(gca, 'ColorOrderIndex', 1);
plot(X1(:,1),X1(:,2),'LineWidth',1.5)
plot(X2(:,1),X2(:,2),'--','LineWidth',1.5)
axis equal;
xlim([-4 4]); ylim([-4 4])
xlabel("$x_1$"); ylabel("$x_2$")
exportgraphics(gca, sprintf('plots/example31_sin2_x_d%i_zoom.pdf',degree), 'ContentType', 'vector');


%% Get system dynamics
degree = 8;
[f, g, h] = getSystem31(degree-1);

%%  Compute the energy functions
eta = 0;
[v] = approxPastEnergy(f, g, h, eta=eta, degree=degree);
[w, K] = approxFutureEnergy(f, g, h, eta=eta, degree=degree);

% Plot the past and future HJB residuals
xPlot = linspace(-1, 1, N); yPlot = linspace(-1, 1, N);
[X, Y] = meshgrid(xPlot, yPlot);

vRES = computeResidualPastHJB(f, g, h, eta, v, degree, 1, N);
wRES = computeResidualFutureHJB(f, g, h, eta, w, degree, 1, N);

figure('Position',[880 287 387 318]);
contourf(X, Y, vRES); colorbar; axis equal
xlabel('$x_1$'); ylabel('$x_2$'); title('Controllability Energy Residual')
zlim([-.8 .6])

figure('Position',[1270 287 387 318]);
contourf(X, Y, wRES); colorbar; axis equal
xlabel('$x_1$'); ylabel('$x_2$'); title('Observability Energy Residual')
zlim([-.3 .4])
drawnow

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
[fbal,gbal,hbal,Tbal] = getBalanceThenReduceRealization(f,g,h,eta=0,transformationDegree=degree-1);
TbalInv = transformationInverse(Tbal);
% FofXU = @(x,u) kronPolyEval(f,x) + kronPolyEval(g,x,scenario='G(x)')*u;
FtofZU = @(z,u) kronPolyEval(fbal,z) + kronPolyEval(gbal,z,scenario='G(x)')*u;

%% Simulate dynamics: u = 5 sin(t/π)
% Solve for z0 initial condition with a Newton type iteration
x0 = [1; 1];
z0 = kronPolyEval(TbalInv,x0);

[t, X1] = ode45(@(t, x) FofXU(x, 1*sin(t/pi)), [0, 40], x0);
[~, Z1] = ode45(@(t, z) FtofZU(z, 1*sin(t/pi)), t, z0);

% Transform Z trajectory into X coordinates to compare
X2 = zeros(size(Z1));
for i = 1:length(Z1)
    X2(i,:) = kronPolyEval(Tbal,Z1(i,:).');
end
Z2 = zeros(size(X1));
for i = 1:length(X1)
    Z2(i,:) = kronPolyEval(TbalInv,X1(i,:).');
end

% Plot grid transformations
% Parameters
numLines = 27; numPoints = 201; lim = 6.75;
% Generate original z coordinates
[yH, xH] = meshgrid(linspace(-lim, lim, numLines), linspace(-2*lim, 2*lim, 2*numPoints)); % Vertical lines
[xV, yV] = meshgrid(linspace(-2*lim, 2*lim, 2*numLines), linspace(-lim, lim, numPoints)); % Horizontal lines

% Compute transformed coordinates
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
zxHtr = zeros(size(xH)); zyHtr = zeros(size(yH));
zxVtr = zeros(size(xV)); zyVtr = zeros(size(yV));
for i=1:length(xH(:))
    [xHtr(i), yHtr(i)] = kronPolyEval(Tbal,[xH(i);yH(i)]);
    [xVtr(i), yVtr(i)] = kronPolyEval(Tbal,[xV(i);yV(i)]);
    
    [zxHtr(i), zyHtr(i)] = kronPolyEval(TbalInv,[xH(i);yH(i)]);
    [zxVtr(i), zyVtr(i)] = kronPolyEval(TbalInv,[xV(i);yV(i)]);
end

% Prepare figure
% Generate figure
figure('Position',[850 442 396 300]); hold on
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i),'Color', [0.75 0.75 0.75])
end
for i=1:size(xV,2)
    plot(xVtr(:,i),yVtr(:,i),'Color', [0.75 0.75 0.75])
end
set(gca, 'ColorOrderIndex', 1);
plot(X1(:,1),X1(:,2),'LineWidth',1.5)
plot(X2(:,1),X2(:,2),'--','LineWidth',1.5)
axis equal;
xlim([-1.75 1.75]); ylim([-1.75 1.75])
xlabel("$x_1$"); ylabel("$x_2$")
exportgraphics(gca, sprintf('plots/example31_sin2_x_d%i.pdf',degree), 'ContentType', 'vector');

figure('Position',[1250 442 396 300]); hold on
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i),'Color', [0.75 0.75 0.75])
end
for i=1:size(xV,2)
    plot(xV(:,i),yV(:,i),'Color', [0.75 0.75 0.75])
end
set(gca, 'ColorOrderIndex', 1);
plot(Z2(:,1),Z2(:,2),'LineWidth',1.5)
plot(Z1(:,1),Z1(:,2),'--','LineWidth',1.5)
axis equal;
xlim([-5 5]); ylim([-5 5])
xlabel("$z_1$"); ylabel("$z_2$")
exportgraphics(gca, sprintf('plots/example31_sin2_z_d%i.pdf',degree), 'ContentType', 'vector');

figure('Position',[290 442 560 300]); hold on
plot(t,X1(:,1),'LineWidth',1.5,'DisplayName','Original realization solution')
plot(t,X2(:,1),'--','LineWidth',1.5,'DisplayName','Balanced realization solution')
legend; grid on
xlabel("Time, t"); ylabel("Angular Position, $x_1$")
exportgraphics(gca, sprintf('plots/example31_sin2_d%i.pdf',degree), 'ContentType', 'vector');



end
