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

fprintf('Running Example 12, polynomial balanced realization...\n')

if nargin < 2
    lim = 3;
    if nargin < 1
        degree = 8;
    end
end

%% Get system dynamics
[f, g, h] = getSystem12(degree - 1, false);  % Scherpen model

%% Compute balanced realization
[fbal,gbal,hbal,Tbal] = getBalanceThenReduceRealization(f,g,h,eta=0,transformationDegree=degree-1);
TbalInv = transformationInverse(Tbal);

fprintf('  - The balanced realization for the nonlinear model is:\n')
dispKronPoly(fbal)
dispKronPoly(gbal)


[v] = approxPastEnergy(f, g, h, eta=0, degree=degree);
[w] = approxFutureEnergy(f, g, h, eta=0, degree=degree);
dispKronPoly(v,n=2),fprintf('\b'),dispKronPoly(w,n=2)

[sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree=degree-1, verbose=true);
[ft,gt,ht] = transformDynamics(f,g,h,TinOd,degree=degree-1);

[vt] = approxPastEnergy(ft, gt, ht, eta=0, degree=degree);
[wt] = approxFutureEnergy(ft, gt, ht, eta=0, degree=degree);
dispKronPoly(vt,n=2),fprintf('\b'),dispKronPoly(wt,n=2)
%% Plot grid transformations
fprintf('   Plotting coordinate grids under the balancing transformation...\n')
% Parameters
numLines = 41; numPoints = 201;

% Generate original coordinates
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

x0 = [1 1].'*(0.2*lim);

% Solve for z0 initial condition with a Newton type iteration
z0 = kronPolyEval(TbalInv,x0);
fprintf(['   Simulating the system in the original vs transformed coordinates for initial condition x0 = [', repmat('%2.2e ', 1, numel(x0)), '], ...\n'], x0)
fprintf(['   ... the transformed initial condition is z0 = [', repmat('%2.2e ', 1, numel(z0)), '] '], z0)
fprintf(' (error: %2.2e). \n', norm(kronPolyEval(Tbal,z0)-x0))

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

fprintf('    -> The figure confirms that the solution trajectories are identical in the original vs transformed coordinates <a href="matlab:nexttile(4), xlim([-3 3])">locally</a>.\n')

end
