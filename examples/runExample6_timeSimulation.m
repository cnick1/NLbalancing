function runExample6_timeSimulation(x0)
%runExample6_timeSimulation Simulate the finite element beam example.
%
%   Usage:  runExample6_timeSimulation(x0)
%
%   Input: x0 - Initial condition scaling factor
%
%   Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” 2023.
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 6\n')

if nargin < 1
    x0 = 0.1;
end

numEls = 3; numNodes = numEls + 1; n = 6*numEls;

% Compute energy functions and CPU time
[f, g, ~, initialCondition] = getSystem6(numEls, 2);
% initialCondition([1,4,7]) = 0;
IC = x0 * initialCondition;
% IC = 0*initialCondition;

tspan  = [0, .01]; 
u = [0; 0]; 
% u = [-5e9; 5e9];
% u = [-1e10; 1e10];
% [t, X] = ode23(@(t, x) [FofXU(f(1), g(1), x(1:n), u); FofXU(f, g, x(n+1:2*n), u)], tspan, [IC; IC]);
[t, X] = ode23(@(t, x) [FofXU(f(1), g(1), x(1:n), u); FofXU(f, g, x(n+1:2*n), u)], tspan, [IC; IC]);

Xlin = X(:,1:n); Xnonlin = X(:,n+1:2*n);

figure
while true
    for i=1:round(length(t)/200):length(t)
        x = linspace(0, 1, numNodes);
        ulin = [0, Xlin(i,1:3:3*numNodes-3)];
        vlin = [0, Xlin(i,2:3:3*numNodes-2)];
        unonlin = [0, Xnonlin(i,1:3:3*numNodes-3)];
        vnonlin = [0, Xnonlin(i,2:3:3*numNodes-2)];
        hold off
        plot(x, 0*vlin, 'k', x, 0*vlin, 'k*');
        hold on;
        plot(x+ulin, vlin, 'r', x+ulin, vlin, 'ro');
        plot(x+unonlin, vnonlin, 'b', x+unonlin, vnonlin, 'b^');
        xlabel('x');
        ylabel('y');
        title(sprintf('%f',t(i)))
        xlim([0 1.25])
        ylim([-.5 .5])
        pause(0.01)
    end
end

end

function xdot = FofXU(f, g, x, u)
% Evaluate the dynamics xdot = f(x, u)
% f = f(1); g = g(1);

FofX = kronPolyEval(f, x);
GofX = g{1}; lg = length(g); Im = speye(size(g{1},2));
xk = 1; for k=2:lg; xk = kron(xk, x); GofX = GofX + g{k} * kron(xk, Im); end

xdot = FofX + GofX*u;
end