function runExample16_balancedReduction(degree)
%runExample16_balancedReduction 
%
%   Usage:  runExample16_balancedReduction(degree,lim)
%
%   Inputs:   
%
%   Description: 
%
%   References: 
%
%   Part of the NLbalancing repository.
%%
% close all;
arguments
    degree = 6
end
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')
fprintf('Running Example 16, polynomial balanced reduced-order model...\n')

%% Get dynamics and define control problem
[f, g, h] = getSystem16(5); n=6;
% g = g(1);
%% Compute balanced realization
r=2;
[fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g,h,r=r,eta=0,degree=5,transformationDegree=degree-1);
TbalInv = transformationInverse(Tbal);

fprintf('  - The balanced realization for the nonlinear model is:\n')
% dispPolyDynamics(fbal,gbal,hbal)
% dispPolyDynamics(fbal2,gbal2,hbal2)

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
G = @(x) kronPolyEval(g, x, scenario='G(x)');
Fbal = @(z) kronPolyEval(fbal, z);
Gbal = @(z) kronPolyEval(gbal, z, scenario='G(x)');
u = @(t) -2e0*sin(.225*t);

%% Simulate original and transformed systems; same input should give same output
x0 = zeros(n,1);
% Solve for z0 initial condition
z0 = kronPolyEval(TbalInv,x0);
% z0 = newtonIteration(x0, @(z) kronPolyEval(Tbal,z), @(z) jcbn(Tbal,z),maxIter=100,verbose=true,z0=kronPolyEval(TbalInv,x0));

fprintf(['   Simulating the system in the original vs transformed coordinates for initial condition x0 = [', repmat('%2.2e ', 1, numel(x0)), '], ...\n'], x0)
fprintf(['   ... the transformed initial condition is z0 = [', repmat('%2.2e ', 1, numel(z0)), '] '], z0)
fprintf(' (error: %2.2e). \n', norm(kronPolyEval(Tbal,z0)-x0))

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
% Simulate both systems
[t1, X1] = ode45(@(t, x) F(x) + G(x)*u(t), [0, 80], x0);
[t2, Z2] = ode45(@(t, z) Fbal(z) + Gbal(z)*u(t), [0, 80], z0);

% Convert z trajectories to x coordinates
X2 = zeros(length(t2),n);
for i=1:length(t2)
    X2(i,:) = kronPolyEval(Tbal,Z2(i,:).');
end

% Convert x trajectories to z coordinates
Z1 = zeros(length(t1),r);
for i=1:length(t1)
    Z1(i,:) = kronPolyEval(TbalInv,X1(i,:).');
end

% Compute the output for the two trajectories
y1 = zeros(1,length(t1));
for i=1:length(t1)
    y1(:,i) = kronPolyEval(h,X1(i,:).');
end
y2 = zeros(1,length(t2));
for i=1:length(t2)
    y2(:,i) = kronPolyEval(hbal,Z2(i,:).');
end
y3 = zeros(1,length(t2));
for i=1:length(t2)
    y3(:,i) = kronPolyEval(h,X2(i,:).');
end

% Plot state trajectories
figure('Position', [827 220 560 420]);
% subplot(3,1,1); hold on;
% plot(t1,X1)
% plot(t2,X2,'--')
% 
% title('Reconstructed state trajectories'); xlabel('Time t')
% ylabel('x_i(t)')
% legend('FOM x_1','FOM x_2','FOM x_3','ROM x_1','ROM x_2','ROM x_3')
% 
% % Plot transformed state trajectories
% subplot(3,1,2); hold on;
% plot(t1,Z1)
% plot(t2,Z2,'--')
% 
% title('Reconstructed state trajectories'); xlabel('Time t')
% ylabel('z_i(t)')
% legend('FOM z_1','FOM z_2','ROM z_1','ROM z_2')

% Plot outputs
% subplot(3,1,3); 
hold on;
plot(t1,y1)
plot(t2,y2,'--')
plot(t2,y3,':')

title('Model output'); xlabel('Time t')
ylabel('y(t)')
legend('FOM output','ROM output using hbar(z)','ROM output h(x)')

fprintf('\n   The output error is: ||yᵣ(t)-y(t)||₂ = %f \n\n', norm(interp1(t2, y2, 0:.1:80) - interp1(t1, y1, 0:.1:80)))

end

