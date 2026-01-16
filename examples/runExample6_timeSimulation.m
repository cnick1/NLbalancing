function runExample6_timeSimulation(x0)
%runExample6_timeSimulation Simulate the finite element beam example.
%
%   Usage:  runExample6_timeSimulation(x0)
%
%   Input: x0 - Initial condition scaling factor
%
%   Reference: [1] N. A. Corbin and B. Kramer, "Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” 2023.
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 6\n')

if nargin < 1
    x0 = 0.01;
end

numEls = 3; numNodes = numEls + 1; n = 6*numEls;

% Compute energy functions and CPU time
[E, f, g, h, initialCondition] = getSystem6(numEls, 1);
g = g(1);
% initialCondition([1,4,7]) = 0;
IC = x0 * initialCondition;
% IC = 0*initialCondition;
% IC2 = IC; IC2([1,4,7]) = 0;
% IC1 = IC; IC1([2,5,8]) = 0;

tspan  = [0, 500];
u = [0; 0]; 
% u = [-5e9; 5e9];
% u = [-1e10; 1e10];
opts = odeset(OutputFcn=@odeprog);
global T0; T0 = tic;
[t, X] = ode23(@(t, x) [FofXU(f(1), g(1), x(1:n), u); FofXU(f, g, x(n+1:2*n), u)], tspan, [IC; IC], opts);
% [t, X] = ode45(@(t, x) [FofXU(f(1), g(1), x(1:n), u); FofXU(f(1), g(1), x(n+1:2*n), u)], tspan, [2*IC1; IC2]);

Xlin = X(:,1:n);
Xnonlin = X(:,n+1:2*n);

figure
plot(abs(fft(Xlin(:,7))))
hold on; plot(abs(fft(Xnonlin(:,8))))

x = linspace(0, 1, numNodes);

% Precompute the node indices for easier indexing later
nodeIndices = 1:3:(3*numNodes-3);

% Create figure
figure
undeformed = animatedline('Color', 'k','LineWidth',1.25); addpoints(undeformed, x, 0*x);
hLin = animatedline('Marker', 'o', 'Color', 'r');
hNonLin = animatedline('Marker', '^', 'Color', 'b');
xlabel('x');
ylabel('y');
xlim([0 1.25]);
ylim([-.05 .05]);
% grid on;
% legend('Linear', 'Linear Marker', 'Nonlinear', 'Nonlinear Marker');

% Set up VideoWriter
v = VideoWriter('plots/animation.avi');
open(v);

% Animation loop
for i = 1:round(length(t)/400):length(t)
    ulin = [0, Xlin(i,nodeIndices)];
    vlin = [0, Xlin(i,nodeIndices + 1)];
    unonlin = [0, Xnonlin(i,nodeIndices)];
    vnonlin = [0, Xnonlin(i,nodeIndices + 1)];
    
    % Clear existing points
    clearpoints(hLin);
    clearpoints(hNonLin);
    
    % Add new points
    addpoints(hLin, x + ulin, vlin);
    addpoints(hNonLin, x + unonlin, vnonlin);
    
    title(sprintf('%f', t(i)));
    drawnow;
    
    % Capture the plot as a frame
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video file
close(v);


end

function xdot = FofXU(f, g, x, u)
% Evaluate the dynamics xdot = f(x, u)
% f = f(1); g = g(1);

FofX = kronPolyEval(f, x);
GofX = g{1}; lg = length(g); Im = speye(size(g{1},2));
xk = 1; for k=2:lg; xk = kron(xk, x); GofX = GofX + g{k} * kron(xk, Im); end

xdot = FofX + GofX*u;
end
