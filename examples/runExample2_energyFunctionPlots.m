function [v, w] = runExample2_energyFunctionPlots(exportPlotData)
%runExample2_energyFunctionPlots Runs the 2D quadratic-bilinear example from [1] to plot energy functions as contour plots
%
%   Usage:  [v,w] = runExample2_energyFunctionPlots(exportPlotData)
%
%   Inputs:
%       exportPlotData   Boolean variable to determine if plots/data are exported
%
%   Outputs:
%       v,w              are coefficients of the past and future energy
%                        function approximations, respectively.
%
%   Reference: [1] Y. Kawano and J. M. A. Scherpen, “Model reduction by
%               differential balancing based on nonlinear hankel operators,”
%               IEEE Transactions on Automatic Control, vol. 62, no. 7,
%               pp. 3293–3308, Jul. 2017, doi: 10.1109/tac.2016.2628201.
%              [2] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗∞
%               energy functions for polynomial control-affine systems,"
%               IEEE Transactions on Automatic Control, pp. 1–13, 2024,
%               doi: 10.1109/tac.2024.3494472
%%
load(fullfile('utils', 'YlGnBuRescaled.mat'))

%% Process inputs
if nargin < 1
    exportPlotData = false;
end
degree = 4;

%% Get model and compute energy functions
[f, g, h] = getSystem2(kawano=true); % use original Kawano model
fprintf('Running Example 2\n')

eta = 0; % values should be between -\infty and 1.

%  Compute the polynomial approximations to the past future energy function
[v] = approxPastEnergy(f, g, h, eta=eta, degree=degree, verbose=true);
[w] = approxFutureEnergy(f, g, h, eta=eta, degree=degree, verbose=true);

%% Plot the past and future energy functions
N = 301;
xPlot = linspace(-1, 1, N); yPlot = linspace(-1, 1, N);
[X, Y] = meshgrid(xPlot, yPlot);

ePast = zeros(N, N); eFuture = zeros(N, N);
for i = 1:N
    for j = 1:N
        x = [X(i, j); Y(i, j)];
        ePast(i, j) = 0.5 * kronPolyEval(v, x, degree);
        eFuture(i, j) = 0.5 * kronPolyEval(w, x, degree);
    end
end
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex','defaultAxesTickLabelInterpreter', 'latex', 'defaulttextinterpreter', 'latex', 'defaultLegendInterpreter', 'latex');

fig1 = figure('Position',[70 287 363.3333 470.6667]);
contourf(X, Y, ePast, 16, 'w'); colormap(flip(YlGnBuRescaled))
xlabel('$x_1$'); ylabel('$x_2$');
xticks(-1:1); yticks(-1:1)
set(gca, 'FontSize', 16)
h = colorbar('FontSize', 16, 'Location','southoutside');
clim([0 80]); set(h, 'ylim', [0 80])


fig2 = figure('Position',[490 287 363.3333 470.6667]);
contourf(X, Y, eFuture, 16, 'w'); colormap(flip(YlGnBuRescaled))
xlabel('$x_1$'); ylabel('$x_2$');
xticks(-1:1); yticks(-1:1)
set(gca, 'FontSize', 16)
h = colorbar('FontSize', 16,'Location','southoutside');
clim([0 1.5]); set(h, 'ylim', [0 1.5])

if exportPlotData
    exportgraphics(fig1, 'plots/example2b.pdf', 'ContentType', 'vector');
    exportgraphics(fig2, 'plots/example2a.pdf', 'ContentType', 'vector');
else
    figure(fig1); title('Past Energy Function')
    figure(fig2); title('Future Energy Function')
end

end
