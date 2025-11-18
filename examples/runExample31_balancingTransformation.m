function runExample31_balancingTransformation(degree,lim)
%runExample31_balancingTransformation Runs the 2D pendulum example to visualize the nonlinear balancing transformations.
%
%   Usage:  runExample31_balancingTransformation(degree,lim)
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
% close all;
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 14\n')

if nargin < 2
    lim = 3.25;
    if nargin < 1
        degree = 4;
    end
end

%% Get system dynamics
[f, g, h, FofXU] = getSystem31(degree-1);

%% Compute balanced realization
[fbal,gbal,hbal,Tbal] = getBalancedRealization(f,g,h,eta=0,transformationDegree=degree-1);
TbalInv = transformationInverse(Tbal);

%% Compute input-normal/output-diagonal realization
[v] = approxPastEnergy(f, g, h, 0, degree);
[w] = approxFutureEnergy(f, g, h, 0, degree);
[~, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree-1);
[finOd,ginOd,hinOd] = transformDynamics(f,g,h,TinOd);
[vbal, wbal] = transformEnergyFunctions(v,w,Tbal);
[vinOd, winOd] = transformEnergyFunctions(v,w,TinOd);

fprintf("\n  - FOM dynamics:\n\n")
dispKronPoly(f)

fprintf("\n  - Balanced dynamics:\n\n")
dispKronPoly(fbal,degree=degree-1)

fprintf("\n  - Energy Functions:\n\n")
dispKronPoly(v,n=2),dispKronPoly(w,n=2)

fprintf("\n  - Input-normal/output-diagonal energy Functions:\n\n")
dispKronPoly(vinOd,n=2),dispKronPoly(winOd,n=2)

fprintf("\n  - Balanced energy Functions:\n\n")
dispKronPoly(vbal,n=2),dispKronPoly(wbal,n=2)

%% Plot grid transformations
% Parameters
numLines = 41; numPoints = 201;

% Generate original coordinates
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-2*lim, 2*lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-2*lim, 2*lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

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

% Generate figure
figure("Position", [185 89 997 830]);
tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile; title(sprintf("x = Φ(z) for z ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
hold on; xlabel("x_1"); ylabel("x_2")
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i))
    plot(xVtr(:,i),yVtr(:,i))
end
axis equal;
nexttile; title(sprintf("z for z ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
hold on; xlabel("z_1"); ylabel("z_2")
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i))
    plot(xV(:,i),yV(:,i))
end
axis equal;
nexttile; title(sprintf("x = Φ(z) for x ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
hold on; xlabel("x_1"); ylabel("x_2")
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i))
    plot(xV(:,i),yV(:,i))
end
axis equal;
nexttile; title(sprintf("z for x ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
hold on; xlabel("z_1"); ylabel("z_2")
for i=1:size(xH,2)
    plot(zxHtr(:,i),zyHtr(:,i))
    plot(zxVtr(:,i),zyVtr(:,i))
end
axis equal;

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
Ft = @(x) kronPolyEval(fbal, x);

x0 = [1 1].'*(0.25*lim);

% Solve for z0 initial condition with a Newton type iteration
z0 = kronPolyEval(TbalInv,x0);
fprintf(['\n         -> Initial condition: z0 = [', repmat('%2.2e ', 1, numel(z0)), '], '], z0)
fprintf('       error: %2.2e \n', norm(kronPolyEval(Tbal,z0)-x0))

% Simulate both systems
[~, X1] = ode45(@(t, x) F(x), [0, 5], x0);
[~, Z] = ode45(@(t, z) Ft(z), [0, 5], z0);

nexttile(1)
plot(X1(:,1),X1(:,2),'g','LineWidth',1.5)
nexttile(3)
plot(X1(:,1),X1(:,2),'g','LineWidth',1.5)
nexttile(2)
plot(Z(:,1),Z(:,2),'r--','LineWidth',1.5)
nexttile(4)
plot(Z(:,1),Z(:,2),'r--','LineWidth',1.5)

%% Transform Z trajectory into X coordinates to compare
X2 = zeros(size(Z));
for i = 1:length(Z)
    X2(i,:) = kronPolyEval(Tbal,Z(i,:).');
end

nexttile(1)
plot(X2(:,1),X2(:,2),'r--','LineWidth',1.5)
nexttile(3)
plot(X2(:,1),X2(:,2),'r--','LineWidth',1.5)
drawnow
end

