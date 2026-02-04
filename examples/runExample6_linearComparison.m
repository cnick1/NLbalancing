function runExample6_linearComparison()
%runExample6_balancedReduction
%
%   Usage:  runExample6_balancedReduction(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description:
%
%   References:
%
%   Part of the NLbalancing repository.
%%
% close all;
degree = 2;
numEls = 1;
U0 = 1e3;
set(groot,'defaultLineLineWidth',2.5,'defaultTextInterpreter','latex')

numNodes = numEls + 1; n = 6*numEls;
fprintf('Running Example 6, n=%d, degree %d...\n',n,degree)

%% Get dynamics and define control problem
[E, f, g, ~, initialCondition] = getSystem6(numEls, 3, false, false);

g = {eye(n)}; h = {eye(n)}; m = size(g{1},2); p = size(h{1},1);
flin = f(1);
F = @(x) kronPolyEval(f, x);
Flin = @(x) kronPolyEval(flin, x);
G = @(x) kronPolyEval(g, x, scenario='G(x)');

% Option 1: Obtain steady-state via time integration
% u = @(t) [zeros(m-2,1); U0*(t>0.0005); 0];
% [t, X1] = ode23s(@(t, x) F(x) + G(x)*u(t), [0 .1], 0*initialCondition);
% y1 = zeros(p,length(t));
% for i=1:length(t)
%     y1(:,i) = kronPolyEval(h, X1(i,:).');
% end
% figure
% plot(y1(n/2-2,:))
%
% x0 = X1(end,:).';


% Option 2: Obtain steady-state Newton iteration for equilibrium point
fsymmetric = f;
for i=2:length(fsymmetric)
    for j = 1:n
        fsymmetric{i}(j,:) = kronMonomialSymmetrize(fsymmetric{i}(j,:),n,i);
    end
end

u = [zeros(m-2,1); U0; 0];
x0 = newtonIteration(-g{1}*u, @(x) kronPolyEval(f, x), @(x) jcbn(fsymmetric, x), maxIter=100);

%% Simulate to compare original dynamics with linearized dynamics
u = @(t) [zeros(m,1)];

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
% Simulate both systems
t = [0:.0001:.05];
[t, X1] = ode23s(@(t, x) F(x) + G(x)*u(t), t, x0);
[t, X2] = ode23s(@(t, x) Flin(x) + G(x)*u(t), t, x0);

%% Compute and plot the outputs
y1 = zeros(p,length(t));
y2 = zeros(p,length(t));
for i=1:length(t)
    y1(:,i) = kronPolyEval(h, X1(i,:).');
    y2(:,i) = kronPolyEval(h, X2(i,:).');
end

figure(10)
close
figure(10)
subplot(2,1,1)
plot(t,y1(n/2-2,:),'DisplayName','FOM output')
hold on
plot(t,y2(n/2-2,:),':','DisplayName','Linearized FOM output')
xlabel('Time, t'); ylabel('$y_1(t)$'); title('Beam tip horizontal displacement ')
legend('Location','southeast'); %ylim(1e-7*[-.35 .05])
grid on

subplot(2,1,2)
plot(t,y1(n/2-1,:),'DisplayName','FOM output')
hold on
plot(t,y2(n/2-1,:),':','DisplayName','Linearized FOM output')
xlabel('Time, t'); ylabel('$y_2(t)$'); title('Beam tip vertical displacement ')
grid on 
set(gcf,"Position",[545 269 774 420])

exportgraphics(gcf, sprintf('plots/example6_linearized.pdf'), 'ContentType', 'vector');

end

