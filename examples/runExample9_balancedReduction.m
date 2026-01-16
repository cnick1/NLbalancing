function runExample9_balancedReduction(degree,lim,reduction)
%runExample9_balancedReduction 
%
%   Usage:  runExample9_balancedReduction(degree,lim)
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
arguments
    degree = 4
    lim = 1
    reduction = true
end
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 9, polynomial balanced reduced-order model...\n')

n = 8;
eps = 0.5;
r = 2;

%% Get dynamics and define control problem
[f, B, Q, ~, y] = getSystem9Neumann(eps, n);
[~, GainsPPR] = ppr(f, B, Q, 1e-1, 2); K = GainsPPR{1};
f{1} = f{1}+B*K;
g = {B}; h = {B.'};

x0 = 0.1 + cos(2*pi*y).*cos(pi*y); % Modified initial condition
x0 = 0.1*x0; 

%% Compute balanced realization
[fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g,h,r=r,eta=0,degree=3,transformationDegree=degree-1);
TbalInv = transformationInverse(Tbal);

fprintf('  - The balanced realization for the nonlinear model is:\n')
dispPolyDynamics(fbal,gbal,hbal)
% dispPolyDynamics(fbal2,gbal2,hbal2)

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
Fbal = @(z) kronPolyEval(fbal, z);

%% Simulate original and transformed systems; same input should give same output
% Solve for z0 initial condition
z0 = kronPolyEval(TbalInv,x0);
% z0 = newtonIteration(x0, @(z) kronPolyEval(Tbal,z), @(z) jcbn(Tbal,z),maxIter=100,verbose=true,z0=kronPolyEval(TbalInv,x0));

fprintf(['   Simulating the system in the original vs transformed coordinates for initial condition x0 = [', repmat('%2.2e ', 1, numel(x0)), '], ...\n'], x0)
fprintf(['   ... the transformed initial condition is z0 = [', repmat('%2.2e ', 1, numel(z0)), '] '], z0)
fprintf(' (error: %2.2e). \n', norm(kronPolyEval(Tbal,z0)-x0))

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
% Simulate both systems
[t1, X1] = ode45(@(t, x) F(x), [0, 5], x0);
[t2, Z2] = ode45(@(t, z) Fbal(z), [0, 5], z0);

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
y1 = zeros(size(t1));
for i=1:length(t1)
    y1(i) = kronPolyEval(h,X1(i,:).');
end
y2 = zeros(size(t2));
for i=1:length(t2)
    y2(i) = kronPolyEval(hbal,Z2(i,:).');
end
y3 = zeros(size(t2));
for i=1:length(t2)
    y3(i) = kronPolyEval(h,X2(i,:).');
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

fprintf('\n   The output error is: ||yᵣ(t)-y(t)||₂ = %f \n\n', norm(interp1(t2, y2, 0:.1:5) - interp1(t1, y1, 0:.1:5)))

%% Manifold figure plot
% dx = 0.05; lim = 2;
% [z1,z2] = meshgrid(-lim:dx:lim,-lim:dx:lim);
% x1 = zeros(size(z1)); x2 = x1; x3 = x1;
% for i=1:numel(z1)
%     [x1(i), x2(i), x3(i)] = kronPolyEval(Tbal,[z1(i);z2(i)]);
% end
% 
% % figure;
% % surf(z1,z2,z3)
% figure;
% surf(x1,x2,x3)
% xlim([-2 3])
% ylim([-2 2])
% zlim([-20 5])
% drawnow
% % return
% 
% hold on;
% plot3(X1(:,1),X1(:,2),X1(:,3),'g',DisplayName='FOM solution')
% plot3(X2(:,1),X2(:,2),X2(:,3),'r',DisplayName='ROM solution')
% 
% if degree == 2
%     fprintf('    -> The figure confirms that the reduced-order solution trajectories on the linear balanced subspace fail to capture the full dynamics. \n\n')
% else
%     fprintf('    -> The figure confirms that the reduced-order solution trajectories on the balanced manifold produce the same output as the original dynamics. \n\n')
% end
%% Simulate the system's response to a random noise input
% Instead of simulating the nonlinear system, it is much faster to simulate
% the linear system and then transform the solution
if false && reduction
    [fl, gl, hl] = getSystem32(transform=false);
    [fl, gl, hl] = getBalanceThenReduceRealization(fl,gl,hl,eta=0,degree=1);
    
    global T0; opts = odeset(OutputFcn=@odeprog); T0 = tic;
    [T, Z] = ode45(@(t, z) (fl{1}*z + gl{1}*randn(1)), [0, 10], [0;0;0],opts);
    X = zeros(size(Z)); Z = 100*Z; % Scale the output to better see the manifold
    for i=1:length(T)
        X(i,:) = [Z(i,1), Z(i,2), Z(i,3) - Z(i,1)^2 - Z(i,2)^2 - Z(i,1)^3]; % Apply the inverse transformation to get the nonlinear solution
    end
    
    hold on;
    plot3(X(:,1),X(:,2),X(:,3),'r')
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3')
    title('Balanced manifold and white noise response')
end
end

