function runExample12_nonlinearBalancing()
%runExample12_balancingTransformation Runs the 2D Fujimoto/Scherpen example
%to compare linear vs nonlinear balancing using the polynomial and Newton
%iteration approaches.
%
%   Usage:  runExample12_balancingTransformation()
%
%   Description: The 2D model from [1-3] is
%     f(x) = -9x₁+6x₁²x₂+6x₂³-x₁⁵-2x₁³x₂²-x₁x₂⁴
%             -9x₂-6x₁³-6x₁x₂²-x₁⁴x₂-2x₁²x₂³-x₂⁵]
%     g(x) = [3√2(9-6x₁x₂+x₁⁴-x₂⁴)/{9+x₁⁴+2x₁²x₂²+x₂⁴},    
%             √2(-9x₁²-27x₂²+6x₁³x₂+6x₁x₂³-(x₁²+x₂²)³)/{9+x₁⁴+2x₁²x₂²+x₂⁴},  
%             √2(27x₁²+9x₂²+6x₁³x₂+6x₁x₂³+(x₁²+x₂²)³/{9+x₁⁴+2x₁²x₂²+x₂⁴},    
%             3√2(9+6x₁x₂-x₁⁴+x₂⁴)/{9+x₁⁴+2x₁²x₂²+x₂⁴}]
%     h(x) = [2√2(3x₁+x₁x₂²+x₂³)(3-x₁⁴-2x₁²x₂²-x₂⁴)/{1+x₁⁴+2x₁²x₂²+x₂⁴};
%             √2(3x₂-x₁³-x₁x₂²)(3-x₁⁴-2x₁²x₂²-x₂⁴)/{1+x₁⁴+2x₁²x₂²+x₂⁴}])
%
%   We compute the energy functions, the input-normal/output-diagonal
%   transformation, and then the true balancing transformation, given by the
%   composition x = ̅Φ(z̄(z̄) = Φ(𝝋(z̄)). We visualize this mapping
%   from the z̄ coordinates to the x coordinates by forming a grid in the
%   z̄ coordinates and mapping that grid to the x coordinates. This is
%   compared for linear and nonlinear transformations, and the polynomial
%   approach is compared with the more direct Newton iteration approach to
%   illustrate the similarity in terms of accuracy but major speedup with
%   the polynomial approach. 
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

%% Get system dynamics
fprintf('  - The nonlinear full-order model is given by:\n');
fprintf(['       f(x) = -9x₁+6x₁²x₂+6x₂³-x₁⁵-2x₁³x₂²-x₁x₂⁴\n' ...
        '               -9x₂-6x₁³-6x₁x₂²-x₁⁴x₂-2x₁²x₂³-x₂⁵]\n' ...
        '       g(x) = [3√2(9-6x₁x₂+x₁⁴-x₂⁴)/{9+x₁⁴+2x₁²x₂²+x₂⁴},    \n' ...
        '               √2(-9x₁²-27x₂²+6x₁³x₂+6x₁x₂³-(x₁²+x₂²)³)/{9+x₁⁴+2x₁²x₂²+x₂⁴},  \n' ...
        '               √2(27x₁²+9x₂²+6x₁³x₂+6x₁x₂³+(x₁²+x₂²)³/{9+x₁⁴+2x₁²x₂²+x₂⁴},    \n' ...
        '               3√2(9+6x₁x₂-x₁⁴+x₂⁴)/{9+x₁⁴+2x₁²x₂²+x₂⁴}]\n' ...
        '       h(x) = [2√2(3x₁+x₁x₂²+x₂³)(3-x₁⁴-2x₁²x₂²-x₂⁴)/{1+x₁⁴+2x₁²x₂²+x₂⁴};\n' ...
        '               √2(3x₂-x₁³-x₁x₂²)(3-x₁⁴-2x₁²x₂²-x₂⁴)/{1+x₁⁴+2x₁²x₂²+x₂⁴}])\n'])   

degree = 5;
[f, g, h] = getSystem12(degree, false);  % Scherpen model

%% Compute balanced realization
[fbal,gbal,hbal,Tbal,SigmaSquared] = getBalanceThenReduceRealization(f,g,h,eta=0,transformationDegree=degree,verbose=true);

%% Plot grid transformations
fprintf('   Plotting coordinate grids under the balancing transformation...\n')
    
if input('Do you want to create the grid plots? 1 or 0 >> ')
    for degree = 1:2:7
        [f, g, h] = getSystem12(degree, false);  % Scherpen model
        [fbal,gbal,hbal,Tbal,SigmaSquared] = getBalanceThenReduceRealization(f,g,h,eta=0,transformationDegree=degree);
        plotGridTransformation(Tbal,transformationInverse(Tbal),degree)
    end

    [v] = approxPastEnergy(f, g, h, degree=degree); [w] = approxFutureEnergy(f, g, h, degree=4);
    [sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree=3, verbose=false);
    plotGridTransformationNewton(TinOd,sigmaSquared)

    for degree = 1:2:5
        [f, g, h] = getSystem12(degree, false);  % Scherpen model
        [v] = approxPastEnergy(f, g, h, degree=degree+1); [w] = approxFutureEnergy(f, g, h, degree=degree+1);
        [sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree=degree, verbose=false);
        [fbal,gbal,hbal,Tbal,SigmaSquared] = getBalanceThenReduceRealization(f,g,h,eta=0,transformationDegree=degree);
        plotGridTransformationComparison(Tbal,transformationInverse(Tbal),degree,TinOd,sigmaSquared)
    end
end
end

function plotGridTransformation(Tbal,TbalInv,degree)

% Parameters
numLines = 25; numPoints = 201; 

% Generate original coordinates
lim = 1.2;
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

% Compute transformed coordinates
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
for i=1:length(xH(:))
    [xHtr(i), yHtr(i)] = kronPolyEval(Tbal,[xH(i);yH(i)]);
    [xVtr(i), yVtr(i)] = kronPolyEval(Tbal,[xV(i);yV(i)]);
end

% Generate figure
figure; hold on; 
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i),'k','LineWidth',1.25)
    plot(xVtr(:,i),yVtr(:,i),'k','LineWidth',1.25)
end
axis equal;
xlim([-.45 .45]); ylim([-.45 .45])
xtemp = xticklabels; ytemp = yticklabels;
xticklabels([]); yticklabels([])
exportgraphics(gca, sprintf('plots/example12_phi_d%i.pdf',degree), 'ContentType', 'vector');
xticklabels(xtemp); yticklabels(ytemp)
title(sprintf("x = Φ(z) for z ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
xlabel("x_1"); ylabel("x_2")

figure; hold on; 
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i),'k','LineWidth',1.25)
    plot(xV(:,i),yV(:,i),'k','LineWidth',1.25)
end
axis equal;
xlim([-.75 .75]); ylim([-.75 .75])
xtemp = xticklabels; ytemp = yticklabels;
xticklabels([]); yticklabels([])
exportgraphics(gca, 'plots/example12_phi.pdf', 'ContentType', 'vector');
xticklabels(xtemp); yticklabels(ytemp)
title(sprintf("z for z ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
xlabel("z_1"); ylabel("z_2")


lim = 0.5;
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

% Compute transformed coordinates
zxHtr = zeros(size(xH)); zyHtr = zeros(size(yH));
zxVtr = zeros(size(xV)); zyVtr = zeros(size(yV));
for i=1:length(xH(:))
    [zxHtr(i), zyHtr(i)] = kronPolyEval(TbalInv,[xH(i);yH(i)]);
    [zxVtr(i), zyVtr(i)] = kronPolyEval(TbalInv,[xV(i);yV(i)]);
end

figure; hold on;
for i=1:size(xH,2)
    plot(xH(:,i),yH(:,i),'k','LineWidth',1.25)
    plot(xV(:,i),yV(:,i),'k','LineWidth',1.25)
end
axis equal;
xlim([-.45 .45]); ylim([-.45 .45])
xtemp = xticklabels; ytemp = yticklabels;
xticklabels([]); yticklabels([])
exportgraphics(gca, 'plots/example12_phiinv.pdf', 'ContentType', 'vector');
xticklabels(xtemp); yticklabels(ytemp)
title(sprintf("x for x ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
xlabel("x_1"); ylabel("x_2")

figure; hold on;
for i=1:size(xH,2)
    plot(zxHtr(:,i),zyHtr(:,i),'k','LineWidth',1.25)
    plot(zxVtr(:,i),zyVtr(:,i),'k','LineWidth',1.25)
end
axis equal;
xlim([-.75 .75]); ylim([-.75 .75])
xtemp = xticklabels; ytemp = yticklabels;
xticklabels([]); yticklabels([])
exportgraphics(gca, sprintf('plots/example12_phiinv_d%i.pdf',degree), 'ContentType', 'vector');
xticklabels(xtemp); yticklabels(ytemp)
title(sprintf("z = Φ⁻¹(z) for x ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
xlabel("z_1"); ylabel("z_2")

end

function plotGridTransformationNewton(TinOd,sigmaSquared)

% Parameters
numLines = 25; numPoints = 201; 

% Generate original coordinates
lim = 1.2;
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

% Compute transformed coordinates
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
for i=1:length(xH(:))
    [xHtr(i), yHtr(i)] = PhiBar([xH(i);yH(i)],TinOd,sigmaSquared);
    [xVtr(i), yVtr(i)] = PhiBar([xV(i);yV(i)],TinOd,sigmaSquared);
end


% Generate figure
figure; hold on;
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i),'k','LineWidth',1.25)
    plot(xVtr(:,i),yVtr(:,i),'k','LineWidth',1.25)
end
axis equal;
xlim([-.45 .45]); ylim([-.45 .45])
xtemp = xticklabels; ytemp = yticklabels;
xticklabels([]); yticklabels([])
exportgraphics(gca, 'plots/example12_phi_newton.pdf', 'ContentType', 'vector');
xticklabels(xtemp); yticklabels(ytemp)
title(sprintf("x = Φ(z) for z ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
xlabel("x_1"); ylabel("x_2")

lim = 0.5;
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

% Compute transformed coordinates
zxHtr = zeros(size(xH)); zyHtr = zeros(size(yH));
zxVtr = zeros(size(xV)); zyVtr = zeros(size(yV));
for i=1:length(xH(:))
    [zxHtr(i), zyHtr(i)] = PhiBarInv2([xH(i);yH(i)],TinOd,sigmaSquared);
    [zxVtr(i), zyVtr(i)] = PhiBarInv2([xV(i);yV(i)],TinOd,sigmaSquared);
end

figure; hold on; 
for i=1:size(xH,2)
    plot(zxHtr(:,i),zyHtr(:,i),'k','LineWidth',1.25)
    plot(zxVtr(:,i),zyVtr(:,i),'k','LineWidth',1.25)
end
axis equal;
xlim([-.75 .75]); ylim([-.75 .75])
xtemp = xticklabels; ytemp = yticklabels;
xticklabels([]); yticklabels([])
exportgraphics(gca, 'plots/example12_phiinv_newton.pdf', 'ContentType', 'vector');
xticklabels(xtemp); yticklabels(ytemp)
title(sprintf("z = Φ⁻¹(z) for x ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
xlabel("z_1"); ylabel("z_2")

end

function plotGridTransformationComparison(Tbal, TbalInv, degree, TinOd,sigmaSquared)

% Parameters
numLines = 25; numPoints = 201; 

% Generate original coordinates
lim = 1.2;
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

% Compute transformed coordinates
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
for i=1:length(xH(:))
    [xHtr(i), yHtr(i)] = kronPolyEval(Tbal,[xH(i);yH(i)]);
    [xVtr(i), yVtr(i)] = kronPolyEval(Tbal,[xV(i);yV(i)]);
end

% Generate figure
figure; hold on;
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i),'k','LineWidth',1.25)
    plot(xVtr(:,i),yVtr(:,i),'k','LineWidth',1.25)
end
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
for i=1:length(xH(:))
    [xHtr(i), yHtr(i)] = PhiBar([xH(i);yH(i)],TinOd,sigmaSquared);
    [xVtr(i), yVtr(i)] = PhiBar([xV(i);yV(i)],TinOd,sigmaSquared);
end
hold on;
for i=1:size(xH,2)
    plot(xHtr(:,i),yHtr(:,i),'r--','LineWidth',1)
    plot(xVtr(:,i),yVtr(:,i),'r--','LineWidth',1)
end


axis equal;
xlim([-.45 .45]); ylim([-.45 .45])
xtemp = xticklabels; ytemp = yticklabels;
xticklabels([]); yticklabels([])
exportgraphics(gca, sprintf('plots/example12_phi_d%i_comparison.pdf',degree), 'ContentType', 'vector');
xticklabels(xtemp); yticklabels(ytemp)
title(sprintf("x = Φ(z) for z ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
xlabel("x_1"); ylabel("x_2")

lim = 0.5;
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

% Compute transformed coordinates
zxHtr = zeros(size(xH)); zyHtr = zeros(size(yH));
zxVtr = zeros(size(xV)); zyVtr = zeros(size(yV));
for i=1:length(xH(:))
    [zxHtr(i), zyHtr(i)] = kronPolyEval(TbalInv,[xH(i);yH(i)]);
    [zxVtr(i), zyVtr(i)] = kronPolyEval(TbalInv,[xV(i);yV(i)]);
end

figure; hold on; 
for i=1:size(xH,2)
    plot(zxHtr(:,i),zyHtr(:,i),'k','LineWidth',1.25)
    plot(zxVtr(:,i),zyVtr(:,i),'k','LineWidth',1.25)
end


for i=1:length(xH(:))
    [zxHtr(i), zyHtr(i)] = PhiBarInv2([xH(i);yH(i)],TinOd,sigmaSquared);
    [zxVtr(i), zyVtr(i)] = PhiBarInv2([xV(i);yV(i)],TinOd,sigmaSquared);
end
hold on; 
for i=1:size(xH,2)
    plot(zxHtr(:,i),zyHtr(:,i),'r--','LineWidth',1)
    plot(zxVtr(:,i),zyVtr(:,i),'r--','LineWidth',1)
end
axis equal;
xlim([-.75 .75]); ylim([-.75 .75])
xtemp = xticklabels; ytemp = yticklabels;
xticklabels([]); yticklabels([])
exportgraphics(gca, sprintf('plots/example12_phiinv_d%i_comparison.pdf',degree), 'ContentType', 'vector');
xticklabels(xtemp); yticklabels(ytemp)
title(sprintf("z = Φ⁻¹(z) for x ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
xlabel("z_1"); ylabel("z_2")

end

function [zbar1, zbar2] = PhiBarInv2(x,TinOd,sigmaSquared)
zbar = newtonIteration(x, @(z) PhiBar(z,TinOd,sigmaSquared), @(z) PhiBarJacobian(z,TinOd,sigmaSquared),maxIter=100);
zbar1 = zbar(1); zbar2 = zbar(2);
end