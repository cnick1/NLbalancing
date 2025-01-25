function runExample24(degree)
%runExample24 Runs the 2D example to visualize the nonlinear transformations.
%
%   Usage:  runExample24(degree)
%
%   Inputs:
%       degree - desired degree of the energy function approximation
%
%   Description: This simple 2D example aims to capture the key idea in the
%   model reduction problem: the presence of a subsystem that in some sense
%   contributes little (perhaps is decoupled) to the overall dynamics, yet
%   drives interactions that cannot directly be eliminated. The model is:
%           xdot_1 = −2 x1 + 20 x1 x2 + u,
%           xdot_2 = −5 x2 + u,
%                y = x1 + x2.
%   Despite being so simple, this is a challenging problem because the
%   nonlinear interaction is strong: those terms are much larger than the
%   linear terms!
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
% close all;
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 24\n')

if nargin < 1
    degree = 4;
end

lim=1;

[f, g, h] = getSystem24(false); eta = 0;
% f{2} = 0.5*f{2}; % For testing scaling f2
% [f, g, h] = getSystem11(5, 1, 10); eta = 1; % Pendulum

F = @(x) kronPolyEval(f, x);

%% Plot open-loop phase portrait in original coordinates
x01 = linspace(-1, 1, 31); x02 = linspace(-1, 1, 31);
figure; hold on; title("Open-loop phase portrait")
for i = 1:length(x01)
    for j = 1:length(x02)
        [t, X] = ode45(@(t, x) F(x), [0, .1], [x01(i);x02(j)]);
        plot(X(:,1),X(:,2))
    end
end
xlim([-1 1]); ylim([-1 1]); xlabel("x_1"); ylabel("x_2");

%%  Compute the energy functions
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~ Computing energy functions:  ~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

% Print energy functions
n = 2; x = sym('x', [1, n]).'; syms(x);

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
fprintf(" ~~~~~~~~~~~ Computing transformation and singular value functions:  ~~~~~~~~~~~~ \n")
[sigmaSquared, Tod] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);

% Print the transformation x = Φ(z)
fprintf("    Transformation: \n        x = \n"); disp(vpa(kronPolyEval(Tod, x), 3))

% Plot the squared singular value functions
z = linspace(- 1, 1, 51);
figure; hold on; title("Singular value functions (squared)")
for i = 1:2
    plot(z, (polyval(flip(sigmaSquared(i, :)), z)))
end
xlabel('z_i'); ylabel('\sigma_i^2'); legend('\sigma_1^2','\sigma_2^2')

%% Plot grid transformations
% Parameters
numLines = 21; numPoints = 201;

% Generate original z coordinates
[xH, yH] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Horizontal lines
[yV, xV] = meshgrid(linspace(-lim, lim, numLines), linspace(-lim, lim, numPoints)); % Vertical lines

% Compute transformed coordinates
xHtr = zeros(size(xH)); yHtr = zeros(size(yH));
xVtr = zeros(size(xV)); yVtr = zeros(size(yV));
for i=1:length(xH(:))
    temp = kronPolyEval(Tod, [xH(i);yH(i)]);
    xHtr(i) = temp(1); yHtr(i) = temp(2);
    temp = kronPolyEval(Tod, [xV(i);yV(i)]);
    xVtr(i) = temp(1); yVtr(i) = temp(2);
end

% Prepare figure
figure("Position", [185 337.6667 997.3333 420]);
subplot(1,2,2); title(sprintf("z for z ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
hold on; axis equal; xlabel("z_1"); ylabel("z_2")
for i=1:numLines
    plot(xH(:,i),yH(:,i))
    plot(xV(:,i),yV(:,i))
end

subplot(1,2,1); title(sprintf("x = Φ(z) for z ∈[-%1.1f,%1.1f] \\times [-%1.1f,%1.1f] ",lim,lim,lim,lim))
hold on; axis equal; xlabel("x_1"); ylabel("x_2")
for i=1:numLines
    plot(xHtr(:,i),yHtr(:,i))
    plot(xVtr(:,i),yVtr(:,i))
end

%% Simulate dynamics
% Compute transformed dynamics
[ft, gt, ht] = transformDynamicsDontInvert(f, g, h, Tod);
[ft2, gt2, ht2] = transformDynamics(f, g, h, Tod);

% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
Ft = @(z) kronPolyEval(ft, z);
Ft2 = @(z) kronPolyEval(ft2, z);

x0 = [1 1].'*(0.1*lim);

% Solve for z0 initial condition with a Newton type iteration
fprintf("    Using Newton iteration to find transformed initial condition ... ")
iter = 0; z0 = Tod{1}\x0; % Initial guess
while iter < 10  % can change max iterations
    iter = iter + 1;
    
    % Update guess iteratively
    z0 = z0 - jcbn(Tod, z0)\(kronPolyEval(Tod, z0)-x0); % Proper Newton iteration
    
    % Break upon convergence within tolerance
    if norm(kronPolyEval(Tod, z0) - x0) < 1e-9 % can edit tolerance
        fprintf("converged in %i iterations. ",iter)
        break;
    end
end
fprintf("\n         -> Initial condition: z0 = [%2.2e %2.2e]. Error: %2.2e \n", z0.', norm(kronPolyEval(Tod, z0) - x0))

% Simulate both systems
[t1, X1] = ode45(@(t, x) F(x), [0, 5], x0);
[t2, Z] = ode45(@(t, z) jcbn(Tod,z)\Ft(z), [0, 5], z0);
[~, Zb] = ode45(@(t, z) Ft2(z), [0, 5], z0);

subplot(1,2,1)
plot(X1(:,1),X1(:,2),'LineWidth',1.5)
subplot(1,2,2)
plot(Z(:,1),Z(:,2),'r','LineWidth',1.5)
plot(Zb(:,1),Zb(:,2),'g--','LineWidth',1.5)


%% Transform Z trajectory into X coordinates to compare
X2 = zeros(size(Z));
for i = 1:length(Z)
    X2(i,:) = kronPolyEval(Tod, Z(i,:).');
end

subplot(1,2,1)
plot(X2(:,1),X2(:,2),'r--','LineWidth',1.5)
end


