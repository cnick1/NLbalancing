function runExample1_regionOfAccuracy_Lyapunov()
%runExample1_regionOfAccuracy_Lyapunov Runs 1D ODE example to compare computed and
%analytical energy functions.
%
%   Usage:  runExample1_regionOfAccuracy_Lyapunov()
%
%   Part of the NLbalancing repository.
%%

% if nargin < 1
%     exportData = false; %change
% end

%% 1st Figure: all energy functions, big mess but just for me.

xd = linspace(-6, 6, 250);
eta = 0.5; % values should be between -\infty and 1.

[f, g, h] = getSystem1();

%  Compute the polynomial approximations to the future energy function
d = 8;
[w] = approxFutureEnergy(f, g, h, eta, d);

w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};

x = linspace(-6, 6, 250);

Ef2 = 0.5 * w2 * x .^ 2;
Ef3 = Ef2 + 0.5 * w3 * x .^ 3;
Ef4 = Ef3 + 0.5 * w4 * x .^ 4;
Ef5 = Ef4 + 0.5 * w5 * x .^ 5;
Ef6 = Ef5 + 0.5 * w6 * x .^ 6;
Ef7 = Ef6 + 0.5 * w7 * x .^ 7;
Ef8 = Ef7 + 0.5 * w8 * x .^ 8;

%   Efd = Ef2;
%   for idx = 3:length(w)
%     Efd = Efd + 0.5*w{idx}*x.^idx;
%   end

%  Compute the analytical solution for comparison

EPlusAnalytic = EgammaPlusNumerical(x, f, g, h, eta);

figure
plot(x(1:10:end), EPlusAnalytic(1:10:end), '+', ...
    x, Ef2, ...
    x, Ef4, ...
    x, Ef6, ...
    x, Ef8, ...
    'LineWidth', 2)
legend('analytic', ...
    'degree 2', ...
    'degree 4', ...
    'degree 6', ...
    'degree 8', ...
    'Location', 'northwest')
xlabel('$x$', ...
    'interpreter', 'latex', ...
    'FontSize', 20, ...
    'fontweight', 'bold')
ylabel('$\mathcal{E}_\gamma^+$', ...
    'interpreter', 'latex', ...
    'FontSize', 20, ...
    'fontweight', 'bold')

ylim([0 20])

regn2 = computeRegion(f, g, h, eta, w, 2, xd, Ef2);
regn4 = computeRegion(f, g, h, eta, w, 4, xd, Ef4);
regn6 = computeRegion(f, g, h, eta, w, 6, xd, Ef6);
regn8 = computeRegion(f, g, h, eta, w, 8, xd, Ef8);

set(gca, 'ColorOrderIndex', 2)
hold on
plot(xd(regn2), regn2(regn2) * 0, '*', 'LineWidth', 1)
plot(xd(regn4), regn4(regn4) * .2, '*', 'LineWidth', 1)
plot(xd(regn6), regn6(regn6) * .4, '*', 'LineWidth', 1)
plot(xd(regn8), regn8(regn8) * .6, '*', 'LineWidth', 1)

tspan = [0, 2]; X0 = -5;
[t2, X2] = ode23(@(t, x) kronPolyEval(f, x) - eta * (g{1} + kronPolyEval(g(2:end), x)) * (g{1} + kronPolyEval(g(2:end), x)).' * (0.5 * kronPolyDerivEval(w(1:2), x)), tspan, X0);
[t4, X4] = ode23(@(t, x) kronPolyEval(f, x) - eta * (g{1} + kronPolyEval(g(2:end), x)) * (g{1} + kronPolyEval(g(2:end), x)).' * (0.5 * kronPolyDerivEval(w(1:4), x)), tspan, X0);
[t6, X6] = ode23(@(t, x) kronPolyEval(f, x) - eta * (g{1} + kronPolyEval(g(2:end), x)) * (g{1} + kronPolyEval(g(2:end), x)).' * (0.5 * kronPolyDerivEval(w(1:6), x)), tspan, X0);
[t8, X8] = ode23(@(t, x) kronPolyEval(f, x) - eta * (g{1} + kronPolyEval(g(2:end), x)) * (g{1} + kronPolyEval(g(2:end), x)).' * (0.5 * kronPolyDerivEval(w(1:8), x)), tspan, X0);

figure; hold on; set(gca, 'ColorOrderIndex', 2)
plot(t2, X2, ...
    t4, X4, ...
    t6, X6, ...
    t8, X8, ...
    'LineWidth', 2)

% figure; hold on;
% for x = xd
%     plot(x, kronPolyEval(f, x),'*')
% end

error("This is wrong! Need Lyapunov sublevelset contained in D! Omega_c ")

end

function regn = computeRegion(f, g, h, eta, w, degree, xd, E)

dataRange = max(xd); N = length(xd);

w = w(1:degree);
n = length(f{1}); % Get state dimension from A matrix

if n > 1 % Remove extra zeros if the model has been tweaked; these should all be matrices
    f = f(~cellfun(@isscalar, f));
    g = g(~cellfun(@isscalar, g));
    h = h(~cellfun(@isscalar, h));
end

% Generate linspace(-dataRange, dataRange, N) once
xn = linspace(-dataRange, dataRange, N);

% Iterate through all points in the discretized space and compute the RES
Lie = zeros(N ^ n, 1);
for i = 1:N ^ n
    % Calculate the indices for each dimension
    indices = mod(floor((i - 1) ./ N .^ (0:(n - 1))), N) + 1;
    x = flip(xn(indices).'); % This is the ith point x in the state-space

    if length(g) > 1
        % Polynomial input
        Lie(i) = (0.5 * kronPolyDerivEval(w, x)) * kronPolyEval(f, x) ...
            - eta * 0.25 * kronPolyDerivEval(w, x) * (g{1} + kronPolyEval(g(2:end), x)) * (g{1} + kronPolyEval(g(2:end), x)).' * kronPolyDerivEval(w, x).';
    else
        % Linear/constant input B
        Lie(i) = (0.5 * kronPolyDerivEval(w, x)) * kronPolyEval(f, x) ...
            - eta * 0.25 * kronPolyDerivEval(w, x) * g{1} * g{1}.' * kronPolyDerivEval(w, x).';

    end
end

if n > 1
    Lie = reshape(Lie, N * ones(1, n));
end

%% compute region where Lie <=0 and E(x) >= 0

regn = (E >= 0 & Lie.' <= 0);

end

function [Ex] = EgammaPlusNumerical(xd, f, g, h, eta)

syms x;

a = -eta / 2 * (g{1} + kronPolyEval(g(2:end), x)) ^ 2;
b = kronPolyEval(f, x);
c = 1/2 * kronPolyEval(h, x) ^ 2;

dEx2 = (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a);
dEx1 = (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a);

% Ex1 = int(dEx1, x);
% Ex2 = int(dEx2, x);

% y = piecewise(x < 0, Ex1, x >= 0, Ex2);

i = 1; Ex = zeros(1, length(xd));
for xi = xd
    if xi < 0
        Ex(i) = vpaintegral(dEx1, x, 0, xi);
    else
        Ex(i) = vpaintegral(dEx2, x, 0, xi);
    end
    i = i + 1;
end

% fplot(y, [-6, 6], '-.+', 'LineWidth', 2)
% hold on
% plot(xd(1:10:end),Ex(1:10:end),'+', 'LineWidth', 1)

end
