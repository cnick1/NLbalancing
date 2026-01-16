function e = runExample17_balancedReduction(degree, n, r, U)
%runExample17_balancedReduction 
%
%   Usage:  runExample17_balancedReduction(degree,lim)
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
    degree = 2
    n = 16
    r = 4
    U = 1e-1
end
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 17, polynomial balanced reduced-order model...\n')


%% Get dynamics and define control problem
[f, B, C] = getSystem17(5, n / 2);
C = C(n/2,:); % B = C.';
g = {B}; h = {C};

%% Compute balanced realization
[fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g,h,r=r,eta=0,degree=5,transformationDegree=degree-1);
TbalInv = transformationInverse(Tbal);

% fprintf('  - The balanced realization for the nonlinear model is:\n')
% dispPolyDynamics(fbal,gbal,hbal)
% dispPolyDynamics(fbal2,gbal2,hbal2)

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
G = @(x) g{1};
Fbal = @(z) kronPolyEval(fbal, z);
Gbal = @(z) gbal{1};% + kronPolyEval(gbal(2:end), z);

%% Simulate original and transformed systems; same input should give same output
x0 = 0*(1:n).'./n; 
% Solve for z0 initial condition
z0 = kronPolyEval(TbalInv,x0);
% z0 = newtonIteration(x0, @(z) kronPolyEval(Tbal,z), @(z) jcbn(Tbal,z),maxIter=100,verbose=true,z0=kronPolyEval(TbalInv,x0));

% fprintf(['   Simulating the system in the original vs transformed coordinates for initial condition x0 = [', repmat('%2.2e ', 1, numel(x0)), '], ...\n'], x0)
% fprintf(['   ... the transformed initial condition is z0 = [', repmat('%2.2e ', 1, numel(z0)), '] '], z0)
% fprintf(' (error: %2.2e). \n', norm(kronPolyEval(Tbal,z0)-x0))

%% Apply reduction by eliminating z3 (set it and its derivative to zero)
% Simulate both systems
tmax = 20; tspan = 0:.01:tmax;

% u = @(t) 1e-1*sin(10*2*pi*t);
u = @(t) [zeros(n-1,1);U*(t>0.000) ];
[t1, X1] = ode45(@(t, x) F(x) + G(x)*u(t), tspan, x0);
[t2, Z2] = ode45(@(t, z) Fbal(z) + Gbal(z)*u(t), tspan, z0);

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

e = norm(interp1(t2, y2, 0:.1:tmax) - interp1(t1, y1, 0:.1:tmax));
fprintf('\n   The output error is: ||yᵣ(t)-y(t)||₂ = %f \n\n', e)

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

end

