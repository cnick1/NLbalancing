function runExample2_regionOfAccuracy_res(exportPlotData)
%runExample2_regionOfAccuracy_res Runs the 2D quadratic-bilinear example from [1] to plot HJB residuals as contour plots
%
%   Usage:  runExample2_regionOfAccuracy_res(exportPlotData)
%
%   Inputs:
%       exportPlotData   Boolean variable to determine if plots/data are exported
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
[f, g, h] = getSystem2(kawano=true);
fprintf('Running Example 2\n')

eta = 0; % values should be between -\infty and 1.

%  Compute the polynomial approximations to the future energy function
[v] = approxPastEnergy(f, g, h, eta=eta, degree=degree, verbose=true);
[w] = approxFutureEnergy(f, g, h, eta=eta, degree=degree, verbose=true);

%% Plot the past and future HJB residuals
N = 301;
xPlot = linspace(-1, 1, N); yPlot = linspace(-1, 1, N);
[X, Y] = meshgrid(xPlot, yPlot);

vRES = computeResidualPastHJB(f, g, h, eta, v, degree, 1, 301);
wRES = computeResidualFutureHJB(f, g, h, eta, w, degree, 1, 301);

set(groot, 'defaultColorbarTickLabelInterpreter', 'latex','defaultAxesTickLabelInterpreter', 'latex', 'defaulttextinterpreter', 'latex', 'defaultLegendInterpreter', 'latex');

fig1 = figure('Position',[880 287 363.3333 470.6667]);
contourf(X, Y, abs(vRES), 16, 'w'); colormap(flip(YlGnBuRescaled))
xlabel('$x_1$'); ylabel('$x_2$');
xticks(-1:1); yticks(-1:1)
h = colorbar('FontSize', 16,'Location','southoutside');
set(gca, 'FontSize', 16)
% clim([0,2.5e-11])
% set(h,'YTick',0:1e-11:2e-11,'TickLabels',{'0','$1\cdot10^{-11}$','$2\cdot10^{-11}$'})
fprintf('The residual of the HJB equation on the unit square is %g\n', norm(vRES, 'inf'));

fig2 = figure('Position',[1270 287 363.3333 470.6667]);
pcolor(X, Y, abs(wRES)); shading interp; colormap(flip(YlGnBuRescaled))
xlabel('$x_1$'); ylabel('$x_2$');
xticks(-1:1); yticks(-1:1)
h = colorbar('FontSize', 16,'Location','southoutside');
set(gca, 'FontSize', 16)
set(h,'YTick',0:2e-16:4e-16,'TickLabels',{'0','$2\cdot10^{-16}$','$4\cdot10^{-16}$'})
fprintf('The residual of the HJB equation on the unit square is %g\n', norm(wRES, 'inf'));

if exportPlotData
    exportgraphics(fig1, 'plots/example2_pastEnergy_residual.pdf', 'ContentType', 'vector');
    exportgraphics(fig2, 'plots/example2_futureEnergy_residual.pdf', 'ContentType', 'vector');
else
    figure(fig1); title('Past HJB Residual')
    figure(fig2); title('Future HJB Residual')
end
end


