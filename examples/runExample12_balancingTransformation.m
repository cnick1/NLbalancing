function runExample12_balancingTransformation(degree,lim)
%runExample12_balancingTransformation Runs the 2D Fujimoto/Scherpen example to visualize the nonlinear balancing transformations.
%
%   Usage:  runExample12_balancingTransformation(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description: The 2D model from [1-3] is
%       f(x) = [-9x₁ + 6x₁² x₂ + 6x₂³ - x₁⁵ - 2x₁³ x₂² - x₁ x₂⁴
%               -9x₂ - 6x₁³ - 6x₁x₂² - x₁⁴ x₂ - 2x₁² x₂³ - x₂⁵]
%       g(x) = [\frac{3√2(9-6x₁x₂+x₁⁴-x₂⁴)}{9+x₁⁴+2x₁²x₂²+x₂⁴},
%                 \frac{√2(-9x₁² - 27 x₂² + 6 x₁³ x₂ + 6 x₁ x₂³ - (x₁² + x₂²)³)}{9+x₁⁴+2x₁²x₂²+x₂⁴};
%               \frac{√2(27x₁²+9x₂²+6x₁³x₂+6x₁x₂³+(x₁²+x₂²)³}{9+x₁⁴+2x₁²x₂²+x₂⁴},
%                 \frac{3√2(9 + 6 x₁ x₂  - x₁⁴ + x₂⁴)}{9+x₁⁴+2x₁²x₂²+x₂⁴}]
%       h(x) = [\frac{2√2(3x₁ + x₁ x₂² + x₂³)(3 - x₁⁴ - 2x₁² x₂² - x₂⁴)}{1 + x₁⁴ + 2 x₁² x₂² + x₂⁴};
%                 \frac{√2(3x₂ - x₁³ - x₁ x₂²)(3 - x₁⁴ - 2 x₁² x₂² - x₂⁴)}{1 + x₁⁴ + 2 x₁² x₂² + x₂⁴}]
%
%   We compute the energy functions, the input-normal/output-diagonal
%   transformation, and then the true balancing transformation, given by the
%   composition x = ̅Φ(z̄(z̄) = Φ(𝝋(z̄)). We visualize this mapping
%   from the z̄ coordinates to the x coordinates by forming a grid in the
%   z̄ coordinates and mapping that grid to the x coordinates.
%
%   References: [1] K. Fujimoto and J. M. A. Scherpen, “Model reduction
%                for nonlinear systems based on the differential
%                eigenstructure of Hankel operators,” in Proceedings of
%                the 40th IEEE Conference on Decision and Control (Cat.
%                No.01CH37228), IEEE, 2001. doi: 10.1109/cdc.2001.980322
%               [2] K. Fujimoto and J. M. A. Scherpen, “Nonlinear
%                input-normal realizations based on the differential
%                eigenstructure of Hankel operators,” IEEE Transactions
%                on Automatic Control, vol. 50, no. 1, pp. 2–18, Jan.
%                2005, doi: 10.1109/tac.2004.840476
%               [3] K. Fujimoto and J. M. A. Scherpen, “Balanced
%                realization and model order reduction for nonlinear
%                systems based on singular value analysis,” SIAM Journal
%                on Control and Optimization, vol. 48, no. 7, pp.
%                4591–4623, Jan. 2010, doi: 10.1137/070695332
%
%   Part of the NLbalancing repository.
%%
% close all;
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 12\n')

if nargin < 2
    lim = 2.5;
    if nargin < 1
        degree = 4;
    end
end

%% Get system dynamics
[f, g, h] = getSystem12(degree - 1, false);  % Scherpen model
f{1} = full(f{1}); g{1} = full(g{1}); h{1} = full(h{1});

%%  Compute the energy functions
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~ Computing energy functions:  ~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
eta = 0;
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
fprintf(" ~~~~~~~~~~~ Computing transformation and singular value functions:  ~~~~~~~~~~~~ \n")
[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);

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
hold on; axis equal; xlabel("x_1"); ylabel("x_2")
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i))
    plot(xVtr(:,i),yVtr(:,i))
end
% xlim([-.2 .45]);ylim([-.2 .45])

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
Ft = @(z) PhiBarJacobian(z,TinOd,sigmaSquared)\kronPolyEval(f, PhiBar(z,TinOd,sigmaSquared));

x0 = [1 1].'*(0.2*lim);

% Solve for z0 initial condition with a Newton type iteration
z0 = newtonIteration(x0, @(z) PhiBar(z,TinOd,sigmaSquared), @(z) PhiBarJacobian(z,TinOd,sigmaSquared),100,true);

% Simulate both systems
[~, X1] = ode45(@(t, x) F(x), [0, 5], x0);
[~, Z] = ode45(@(t, z) Ft(z), [0, 5], z0);

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
end


