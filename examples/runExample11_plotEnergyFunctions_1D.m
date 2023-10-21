function [v, w] = runExample11_plotEnergyFunctions_1D(exportPlotData, nGterms, degree, eta, varargin)
%runExample5 Runs the 3D unicycle example
%
%   Usage:  [v,w] = runExample11(degree,plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel)
%
%   runExample11() runs the default case of a quadratic model from [1] which
%                 is based on a model from [2].
%
%   Inputs:
%       degree          is the degree of energy function approximations
%       exportPlotData   Boolean variable to determine if plots/data are exported
%
%
%   Outputs:
%       v,w              are coefficients of the past and future energy
%                        function approximations, respectively.
%
%   The value of eta is set below.
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%% Process inputs
degree = 8;

% Compute energy functions
eta = 1; % values should be between -\infty and 1.
% eta=1 is HJB/closed-loop balancing, 0 is open loop.

%% Get model and compute energy functions
m = 1; L = 10; %56.5962*scale;
[f, g, h] = getSystem11(9, m, L);
[w] = pqr(f, g, h, 1 / eta, degree, true);

%% 1D plots
xrange = pi; nPoints = 100;

figure;

% Iterate through each axis and plot
n = length(f{1});
for xi = 1:n
    xs = zeros(nPoints, n);
    xs(:, xi) = linspace(-xrange, xrange, nPoints).';

    energies = zeros(nPoints, 1);
    for i = 1:nPoints
        x = xs(i, :).';
        energies(i) = 0.5 * kronPolyEval(w, x, degree);
    end

    subplot(1, n, xi)
    plot(xs(:, 1), energies, 'LineWidth', 2)

    xlabel(sprintf('$x_%i$', xi), 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold'); ylabel('$\mathcal{E}_\gamma^+$', 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold')

end

%% 2D plots
nX = 301; nY = nX;
xLim = pi; yLim = 5;
xPlot = linspace(-xLim, xLim, nX);
yPlot = linspace(-yLim, yLim, nY);
eFuture = zeros(nY, nX);
[X, Y] = meshgrid(xPlot, yPlot);

for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        eFuture(i, j) = 0.5 * kronPolyEval(w, x, degree);
        if eFuture(i, j) < 0
            eFuture(i, j) = NaN;
        end
    end
end

%     set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
%     set(groot, 'defaulttextinterpreter', 'latex');
%     set(groot, 'defaultLegendInterpreter', 'latex');

fig1 = figure;
% ('Position',[600 50 1000 400])
% subplot(1,2,1)
contourf(X, Y, eFuture, 16, 'w'); hold on;
% logMaxEFuture = log10(max(max(eFuture)));
% contour(X, Y, eFuture, [0, logspace(-3, ceil(logMaxEFuture), 20)] ./ (10 ^ (ceil(logMaxEFuture) - logMaxEFuture)))
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 16)
xticks([-pi, 0, pi])
xticklabels({'-\pi', '0', '\pi'})
%     axis equal

%         set(h, 'ylim', [0 1.5])
load(fullfile('utils', 'YlGnBuRescaled.mat'))
colormap(flip(YlGnBuRescaled))
