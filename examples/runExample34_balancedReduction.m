function runExample34_balancedReduction(degree)
%runExample34_balancedReduction Runs the 2D academic example 
%
%   Part of the NLbalancing repository.
%%

set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 34, polynomial balanced realization...\n')



%% Get system dynamics
[f, g, h] = getSystem34();
x0 = [.1;.1;.1]*0;
%% Compute balanced realization
[fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g,h,r=2,eta=0,degree=3,transformationDegree=degree-1,verbose=true);
TbalInv = transformationInverse(Tbal);

fprintf('  - The balanced realization for the nonlinear model is:\n')
dispPolyDynamics(fbal,gbal,hbal,variable='z')
% dispPolyDynamics(fbal2,gbal2,hbal2)

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
G = @(x) kronPolyEval(g, x, scenario='G(x)');
Fbal = @(z) kronPolyEval(fbal, z);
Gbal = @(z) kronPolyEval(gbal, z, scenario='G(x)');
u = @(t) [1*sin(2*pi*t/1); 0];

%% Simulate original and transformed systems; same input should give same output
% Solve for z0 initial condition
z0 = kronPolyEval(TbalInv,x0);
fprintf(['   Simulating the system in the original vs transformed coordinates for initial condition x0 = [', repmat('%2.2e ', 1, numel(x0)), '], ...\n'], x0)
fprintf(['   ... the transformed initial condition is z0 = [', repmat('%2.2e ', 1, numel(z0)), '] '], z0)
fprintf(' (error: %2.2e). \n', norm(kronPolyEval(Tbal,z0)-x0))

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
% Simulate both systems
[t1, X1] = ode45(@(t, x) F(x) + G(x)*u(t), [0, 5], x0);
[t2, Z2] = ode45(@(t, z) Fbal(z) + Gbal(z)*u(t), [0, 5], z0);

% Compute the output for the two trajectories
p = size(h{1},1);
y1 = zeros(length(t1),p);
y2 = zeros(length(t2),p);
for i=1:length(t1) 
    y1(i,:) = kronPolyEval(h, X1(i,:).');
end
for i=1:length(t2) 
    y2(i,:) = kronPolyEval(hbal, Z2(i,:).');
end

% Plot state trajectories
figure('Position', [827 220 560 420]);
hold on;
plot(t1,y1(:,2))
plot(t2,y2(:,2),'--')

title('Model output'); xlabel('Time t')
ylabel('y(t)')
legend('FOM output','ROM output using hbar(z)','ROM output h(x)')

fprintf('\n   The output error is: ||yᵣ(t)-y(t)||₂ = %f \n\n', norm(interp1(t2, y2(:,2), 0:.1:5) - interp1(t1, y1(:,2), 0:.1:5)))


end