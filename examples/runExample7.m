function [v, w] = runExample7()
%runExample7 Runs 3D aircraft stall model
%   Usage:  [v, w] = runExample7()
%
%   runExample7() runs the aircraft stall stabilization model from
%   Garrard 1977 [1].
%
%   Outputs:
%       v,w             are coefficients of the past and future energy
%                       function approximations, respectively.
%
%   Reference: [1] W. L. Garrard and J. M. Jordan, “Design of nonlinear
%               automatic flight control systems,” Automatica, vol. 13,
%               no. 5, pp. 497–505, Sep. 1977,
%               doi: 10.1016/0005-1098(77)90070-x
%
%   Part of the NLbalancing repository.
%%

[f, g, h] = getSystem7();
% g = g(1);
fprintf('Running Example 7\n')

eta = 1; % values should be between -\infty and 1.
fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

%  Compute the polynomial approximations to the future energy function
d = 8;
% [v] = approxPastEnergy(f, f{2}, g, h, eta, d);
[w] = approxFutureEnergy(f, f{2}, g, h, eta, d);

% syms x1 x2 x3
% vpa(-eta*g{1}.'*(kronPolyDerivEval(w(1:2), [x1;x2;x3]).'/2),3)
% vpa(-eta*g{1}.'*(kronPolyDerivEval(w(1:4), [x1;x2;x3]).'/2),3)

%% Define dynamics and controllers
F = @(x) kronPolyEval(f, x);
G = @(x) (g{1} + kronPolyEval(g(2:end), x));
U2 = @(x) [.47 * x(1); 0; 46]; U3 = [.63; 0; 61.4];
% U2 = 0; U3 = 0;
U2 = @(x) [.47 * x(1); 0; 46]; U3 = [.63; 0; 61.4];

uLinGarrard = @(x) (- 0.053 * x(1) + 0.5 * x(2) + 0.521 * x(3));
uQuadGarrard = @(x) (- 0.053 * x(1) + 0.5 * x(2) + 0.521 * x(3) + 0.04 * x(1) ^ 2 - 0.048 * x(1) * x(2));
uCubGarrard = @(x) (- 0.053 * x(1) + 0.5 * x(2) + 0.521 * x(3) + 0.04 * x(1) ^ 2 - 0.048 * x(1) * x(2) + 0.374 * x(1) ^ 3 - 0.312 * x(1) ^ 2 * x(2));
uLin = @(x) (- eta * G(x).' * kronPolyDerivEval(w(1:2), x).' / 2);
uQuad = @(x) (- eta * G(x).' * kronPolyDerivEval(w(1:3), x).' / 2);
uCub = @(x) (- eta * G(x).' * kronPolyDerivEval(w(1:4), x).' / 2);

%% Recreate Garrard Figure 1
tspan = [0, 5];

% First IC
alpha0 = 22.9;
X0 = [0.0174533 * alpha0; 0; 0];

u = uLinGarrard; [t1, X1] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);
u = uQuadGarrard; [t2, X2] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);
u = uCubGarrard; [t3, X3] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);

figure; hold on;
plot(t1, X1(:, 1) / 0.0174533)
plot(t2, X2(:, 1) / 0.0174533)
plot(t3, X3(:, 1) / 0.0174533)

set(gcf, 'Position', 1.0e+03 * [2.5150 -0.0310 1.0093 0.5927])

% Second IC
alpha0 = 26;
% alpha0 = 25.628;
X0 = [0.0174533 * alpha0; 0; 0];

u = uLinGarrard; [t1, X1] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);
u = uQuadGarrard; [t2, X2] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);
u = uCubGarrard; [t3, X3] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);

set(gca, 'ColorOrderIndex', 1)
hold on;
plot(t1, X1(:, 1) / 0.0174533)
plot(t2, X2(:, 1) / 0.0174533)
plot(t3, X3(:, 1) / 0.0174533)

ylim([0 35])
xlim(tspan)

%%
% set(gca,'ColorOrderIndex',1)
%
% plot(t4,X4(:,1)/0.0174533,':','LineWidth', 2)
% plot(t5,X5(:,1)/0.0174533,':','LineWidth', 2)
% plot(t6,X6(:,1)/0.0174533,':','LineWidth', 2)
% ylim([0 alpha0 + 3])
% xlim(tspan)

end
