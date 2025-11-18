function runExample13_balancingTransformation(degree,lim)
%runExample13_balancingTransformation Runs the 2D example from Gray & Scherpen 2001 [1] to visualize the nonlinear balancing transformations.
%
%   Usage:  runExample13_balancingTransformation(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description: The 2D model from [1] is
%           f(x) = -[α² x₁ + 2 α x₂ + (α² - 2)x₂²;
%                                x₂              ]
%           g(x) = √2[α - 2 x₂;
%                        1    ]
%           h(x) = 1/√3 [3 α (x₁ + x₂²) + (α - 2√2)x₂]
%
%   where α = (√3 + √2)(√3 + 2). We compute the energy functions, the
%   input-normal/output-diagonal transformation, and then the true balancing
%   transformation, given by the composition x = ̅Φ(z̄(z̄) = Φ(𝝋(z̄)). We visualize
%   this mapping from the z̄ coordinates to the x coordinates by forming a grid
%   in the z̄ coordinates and mapping that grid to the x coordinates.
%
%   References: [1] W. S. Gray and J. M. A. Scherpen, “On the nonuniqueness
%               of singular value functions and balanced nonlinear
%               realizations,” Systems & Control Letters, vol. 44, no. 3,
%               pp. 219–232, Oct. 2001, doi: 10.1016/s0167-6911(01)00144-x
%
%   Part of the NLbalancing repository.
%%
% close all;
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 13, polynomial balanced realization...\n')

if nargin < 2
    lim = 1;
    if nargin < 1
        degree = 4;
    end
end

fprintf("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
fprintf(" ~~~~~~~~~~~~~~~~ Beginning comparisons with Gray/Scherpen 2001: ~~~~~~~~~~~~~~~~ \n")
%% Get system dynamics
[f, g, h] = getSystem13();

fprintf('  - The dynamics in the original coordinates are:\n')
dispPolyDynamics(f,g,h)

if degree == 4
fprintf("  - Comparing our energy function with Gray/Scherpen 2001 Example 2.1:\n")
[v] = approxPastEnergy(f, g, h, 0, degree);
[w] = approxFutureEnergy(f, g, h, 0, degree);
fprintf("  - Energy Functions:\n")
dispKronPoly(v,n=2),fprintf('\b'),dispKronPoly(w,n=2)
fprintf("                             ->  Energy functions match.\n\n")
end
%% Compare the first transformation
% The first transformation given is 
%   x = 𝜙(z) = [z₁ + z₁²; z₂]
% It should be
%   x = 𝜙(z) = [z₁ - z₂²; z₂]
fprintf(['  - Validating the first transformation given in Gray/Scherpen 2001 Example 2.1: The \n' ...
    '    first transformation given is x = 𝜙(z) = [z₁ + z₁²; z₂]; it appears it should\n' ...
    '    instead be x = 𝜙(z) = [z₁ - z₂²; z₂]. Furthermore, this is actually the 𝘪𝘯𝘷𝘦𝘳𝘴𝘦\n' ...
    '    transformation, so it should be written as z = 𝜙⁻¹(x) = [x₁ - x₂²; x₂]\n'])
x = sym('x', [1, 2]).'; 
[Tnl,~,~] = approxPolynomialDynamics([x(1) - x(2)^2; x(2)], [1;1], x(1), x, 2);
[ftr,gtr,htr] = transformDynamics(f,g,h,Tnl,degree=5);

fprintf('  - The dynamics in these coordinates are:\n')
dispPolyDynamics(ftr,gtr,htr)

[v] = approxPastEnergy(ftr, gtr, htr, 0, degree);
[w] = approxFutureEnergy(ftr, gtr, htr, 0, degree);

fprintf("  - Energy Functions:\n        ")
dispKronPoly(v,n=2),fprintf('\b        '),dispKronPoly(w,n=2)

fprintf("                             ->  Energy functions match.\n")
fprintf("\n       ->  𝘐𝘵 𝘢𝘱𝘱𝘦𝘢𝘳𝘴 𝘵𝘩𝘢𝘵 𝘵𝘩𝘪𝘴 𝘮𝘰𝘥𝘦𝘭 𝘸𝘢𝘴 𝘤𝘰𝘯𝘴𝘵𝘳𝘶𝘤𝘵𝘦𝘥 𝘧𝘳𝘰𝘮 𝘢\n           𝘭𝘪𝘯𝘦𝘢𝘳 𝘮𝘰𝘥𝘦𝘭 𝘧𝘰𝘭𝘭𝘰𝘸𝘦𝘥 𝘣𝘺 𝘢 𝘯𝘰𝘯𝘭𝘪𝘯𝘦𝘢𝘳 𝘵𝘳𝘢𝘯𝘴𝘧𝘰𝘳𝘮𝘢𝘵𝘪𝘰𝘯! \n\n")




fprintf('  - The balanced realization for this linear model is:\n')
[fbal1,gbal1,hbal1,~] = getBalancedRealization(ftr,gtr,htr,eta=0,degree=1);
dispPolyDynamics(fbal1,gbal1,hbal1)


fprintf(" ~~~~~~~~~~~~~~~~ Finished comparisons with Gray/Scherpen 2001: ~~~~~~~~~~~~~~~~~ \n")
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

%% Compute balanced realization
[fbal,gbal,hbal,Tbal] = getBalancedRealization(f,g,h,eta=0,transformationDegree=degree-1);
TbalInv = transformationInverse(Tbal);

fprintf('  - The balanced realization for the nonlinear model is:\n')
dispPolyDynamics(fbal,gbal,hbal)

if degree == 4
    fprintf(['\n       ->  𝘚𝘪𝘯𝘤𝘦 𝘵𝘩𝘦 𝘣𝘢𝘭𝘢𝘯𝘤𝘦𝘥 𝘳𝘦𝘢𝘭𝘪𝘻𝘢𝘵𝘪𝘰𝘯 𝘪𝘴 𝘭𝘪𝘯𝘦𝘢𝘳 𝘧𝘰𝘳 𝘵𝘩𝘪𝘴 𝘮𝘰𝘥𝘦𝘭,\n ' ...
        '          𝘵𝘩𝘦 𝘣𝘢𝘭𝘢𝘯𝘤𝘪𝘯𝘨 𝘵𝘳𝘢𝘯𝘴𝘧𝘰𝘳𝘮𝘢𝘵𝘪𝘰𝘯 𝘢𝘤𝘵𝘶𝘢𝘭𝘭𝘺 𝘧𝘪𝘯𝘥𝘴 𝘵𝘩𝘦 𝘯𝘰𝘯𝘭𝘪𝘯𝘦𝘢𝘳\n ' ...
        '          𝘪𝘯𝘷𝘦𝘳𝘴𝘦 𝘵𝘳𝘢𝘯𝘴𝘧𝘰𝘳𝘮𝘢𝘵𝘪𝘰𝘯 𝘵𝘰 𝘮𝘢𝘬𝘦 𝘵𝘩𝘦 𝘥𝘺𝘯𝘢𝘮𝘪𝘤𝘴 𝘭𝘪𝘯𝘦𝘢𝘳.\n\n'])
else
    fprintf('\n       ->  The linear transformation fails to find the true balanced realization.\n\n')
end
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

x0 = [1 1].'*(0.5*lim);

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

fprintf('    -> The figure confirms that the solution trajectories are identical in the original vs transformed coordinates.\n')

end


