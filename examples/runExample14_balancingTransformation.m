function runExample14_balancingTransformation(degree,lim)
%runExample14_balancingTransformation Runs the 2D gradient double pendulum [1-3] example to visualize the nonlinear balancing transformations.
%
%   Usage:  runExample14_balancingTransformation(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description: This model has been used several times in the literature [1-3].
%   The system describes a set of dynamics related to the double pendulum;
%   however, where the double pendulum would have 4D dynamics and only
%   marginal stability, the associated 2D gradient system is asymptotically
%   stable, and hence the model was more approachable when method were more
%   limited. The gradient system dynamics are
%       ẋ = -M⁻¹(x) 𝜕V(x)/𝜕x + M⁻¹(x)[1;0] u
%       y = x₁
%   The mass matrix and its inverse are then
%       M(x) = [m₁₁, m₁₂;    M⁻¹(x) = _______1̲_______   [m₂₂, -m₂₁;
%               m₂₁, m₂₂]            (m₁₁m₂₂ - m₁₂m₂₁)  -m₁₂,  m₁₁]
%   where the entries are
%       m₁₁       = m₁ l₁² + m₂ l₁² + m₂ l₂² + 2 m₂ l₁ l₂ cos x₂
%       m₁₂ = m₂₁ = m₂ l₂² + m₂ l₁ l₂ cos x₂
%       m₂₂       = m₂ l₂²
%   The potential energy of the system is
%       V(x) = - m₁ g l₁ cos x₁ - m₂ g (l₁ cos x₁ + l₂ cos(x₁ + x₂))
%
%   We compute the energy functions, the input-normal/output-diagonal
%   transformation, and then the true balancing transformation, given by the
%   composition x = Φbar(z̄) = Φ(𝝋(z̄)). We visualize this mapping
%   from the z̄ coordinates to the x coordinates by forming a grid in the
%   z̄ coordinates and mapping that grid to the x coordinates.
%
%   References: [1] J. M. A. Scherpen, “Balancing for nonlinear systems,”
%               PhD Dissertation, University of Twente, 1994.
%               [2] W. S. Gray and J. M. A. Scherpen, “Minimality and local
%               state decompositions of a nonlinear state space realization
%               using energy functions,” IEEE Transactions on Automatic
%               Control, vol. 45, no. 11, pp. 2079–2086, 2000, doi:
%               10.1109/9.887630
%               [3] K. Fujimoto and J. M. A. Scherpen, “Nonlinear
%               input-normal realizations based on the differential
%               eigenstructure of Hankel operators,” IEEE Transactions on
%               Automatic Control, vol. 50, no. 1, pp. 2–18, Jan. 2005,
%               doi: 10.1109/tac.2004.840476
%
%   Part of the NLbalancing repository.
%%
% close all;
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 14\n')

if nargin < 2
    lim = .5;
    if nargin < 1
        degree = 4;
    end
end

%% Get system dynamics
[f, g, h] = getSystem14(degree - 1, 1);

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

x0 = [1 1].'*(0.5*lim);

% Solve for z0 initial condition with a Newton type iteration
z0 = newtonIteration(x0, @(z) PhiBar(z,TinOd,sigmaSquared), @(z) PhiBarJacobian(z,TinOd,sigmaSquared));


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


