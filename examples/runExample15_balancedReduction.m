function runExample15_balancedReduction(TransformationDegree, DynamicsDegree, r, U0)
%runExample15_balancedReduction 
%
%   Usage:  runExample15_balancedReduction(degree,r)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                   r - reduced order dimension
%
%   Description: This model has been used several times in the literature [1,2].
%   The 4D state-space model is
%       ẋ = [                  x₃;
%                              x₄;
%              M⁻¹( ∂L/∂(x₁x₂) - Ṁ [x₃; x₄] + Q) ]
%       y = [l₁ * sin(x₁) + l₂ * sin(x₁ + x₂);
%            l₁ * cos(x₁) + l₂ * cos(x₁ + x₂)]
%   The output is taken as the horizontal and vertical positions of the
%   second mass. The mass matrix and its inverse are
%       M(x) = [m₁₁, m₁₂;    M⁻¹(x) = _______1_______  [m₂₂, -m₂₁;
%               m₂₁, m₂₂]            (m₁₁m₂₂ - m₁₂m₂₁) -m₁₂,  m₁₁]
%    where the entries are
%       m₁₁       = m₁ l₁² + m₂ l₁² + m₂ l₂² + 2 m₂ l₁ l₂ cos x₂
%       m₁₂ = m₂₁ = m₂ l₂² + m₂ l₁ l₂ cos x₂
%       m₂₂       = m₂ l₂²
%   The Lagrangian is L(x,ẋ) = T(ẋ) - V(x), where the potential energy is
%       V(x) = - m₁ g l₁ cos x₁ - m₂ g (l₁ cos x₁ + l₂ cos(x₁ + x₂))
%   and the kinetic energy is T(ẋ) = 1/2 ẋ.' * M * ẋ.
%
%   We compute the energy functions, the input-normal/output-diagonal
%   transformation, and then the true balancing transformation, given by the
%   composition x = ̅Φ(z̄(z̄) = Φ(𝝋(z̄)). We then simulate the transformed and
%   optionally reduced system and compare with the original (full-order) model.
%
%   References: [1]  K. Fujimoto and D. Tsubakino, “On computation of
%               nonlinear balanced realization and model reduction,” in
%               2006 American Control Conference, IEEE, 2006. doi:
%               10.1109/acc.2006.1655399
%               [2] K. Fujimoto and D. Tsubakino, “Computation of nonlinear
%               balanced realization and model reduction based on Taylor
%               series expansion,” Systems & Control Letters, vol. 57, no.
%               4, pp. 283–289, Apr. 2008, doi: 10.1016/j.sysconle.2007.08.015
%
%   Part of the NLbalancing repository.
%%
% close all;
arguments
    TransformationDegree = 3
    DynamicsDegree = 3
    r = 2
    U0 = 1
end

set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 15, polynomial balanced reduced-order model...\n')
n = 4;

%% Get dynamics and define control problem
[f, g, h, fsym, gsym, hsym] = getSystem15(DynamicsDegree); [p,n] = size(h{1});
fprintf('  - The full-order nonlinear model is:\n')
dispPolyDynamics(f,g,h)

%% Compute balanced realization
[fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g,h,eta=0,degree=DynamicsDegree,transformationDegree=TransformationDegree);
fprintf('  - The full-order balanced realization for the nonlinear model is:\n')
dispPolyDynamics(fbal,gbal,hbal,variable='z')
% plotSingularValueFunctions(sigmaSquared)


[fbal,gbal,hbal,Tbal,sigmaSquared] = getBalanceThenReduceRealization(f,g,h,r=r,eta=0,degree=DynamicsDegree,transformationDegree=TransformationDegree); 
TbalInv = transformationInverse(Tbal);
fprintf('  - The reduced-order balanced realization for the nonlinear model is:\n')
dispPolyDynamics(fbal,gbal,hbal,variable='z')

%% Simulate dynamics
% Compare original dynamics with transformed dynamics
F = @(x) kronPolyEval(f, x);
G = @(x) kronPolyEval(g, x, scenario='G(x)');
Fbal = @(z) kronPolyEval(fbal, z);
Gbal = @(z) kronPolyEval(gbal, z, scenario='G(x)');
m = size(g{1},2);
u = @(t) [zeros(m-1,1); U0*sin(2.5*t)];

%% Simulate original and transformed systems; same input should give same output
% Solve for z0 initial condition
x0 = [.1 .1 .1 .1].'*0;
z0 = kronPolyEval(TbalInv,x0);
% z0 = newtonIteration(x0, @(z) kronPolyEval(Tbal,z), @(z) jcbn(Tbal,z),maxIter=100,verbose=true,z0=kronPolyEval(TbalInv,x0));

fprintf(['   Simulating the system in the original vs transformed coordinates for initial condition x0 = [', repmat('%2.2e ', 1, numel(x0)), '], ...\n'], x0)


%% Apply reduction by eliminating z3 (set it and its derivative to zero)
% Simulate both systems
global T0; 
opts = odeset(OutputFcn=@odeprog);
T0 = tic; [t1, X1] = ode45(@(t, x) F(x) + G(x)*u(t), [0:.1:20], x0, opts);


fprintf(['   ... the transformed initial condition is z0 = [', repmat('%2.2e ', 1, numel(z0)), '] '], z0)
fprintf(' (error: %2.2e). \n', norm(kronPolyEval(Tbal,z0)-x0))
clear opts;
opts = odeset(OutputFcn=@odeprog);
T0 = tic; [t2, Z2] = ode45(@(t, z) Fbal(z) + Gbal(z)*u(t), [0:.1:20], z0, opts);

%% Compute and plot the outputs 
p = size(h{1},1);
y1 = zeros(p,length(t1));
y2 = zeros(p,length(t2));
for i=1:length(t1) 
    y1(:,i) = kronPolyEval(h, X1(i,:).');
end
for i=1:length(t2) 
    y2(:,i) = kronPolyEval(hbal, Z2(i,:).');
end

figure
plot(y1(1,:))
hold on
plot(y2(1,:),':')

figure
plot(y1(2,:))
hold on
plot(y2(2,:),':')

% Plot state trajectories
% figure('Position', [827 220 560 420]);
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
% hold on;
% plot(t1,y1)
% plot(t2,y2,'--')
% plot(t2,y3,':')
% 
% title('Model output'); xlabel('Time t')
% ylabel('y(t)')
% legend('FOM output','ROM output using hbar(z)','ROM output h(x)')
% 
% fprintf('\n   The output error is: ||yᵣ(t)-y(t)||₂ = %f \n\n', norm(interp1(t2, y2, 0:.1:5) - interp1(t1, y1, 0:.1:5)))

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

end

