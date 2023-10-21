function [v, w] = runExample11_plotEnergyFunctions(exportPlotData, nFterms, degree, eta, varargin)
%runExample11_plotEnergyFunctions Runs the 2D example to plot energy functions as contour plots
%
%   Usage:  [v,w] = runExample11_plotEnergyFunctions(degree,plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel)
%
%   runExample11_plotEnergyFunctions() runs the default case of a quadratic model from [1] which
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
if nargin < 4
    if nargin < 3
        if nargin < 2
            if nargin < 1
                exportPlotData = false;
            end
            nFterms = 1;
        end
        degree = nFterms + 1;
    end
    % Compute energy functions
    eta = 1; % values should be between -\infty and 1.
    % eta=1 is HJB/closed-loop balancing, 0 is open loop.
end

if nFterms == 1
    nFterms = 2; % Note F2 is zero; this is just to be able to compute a controller and ignore the error if F2 doesn't exist
end

%% Get model and compute energy functions
scale = .1767; scaling = 1 / sqrt(scale); % For plot and initial condition scaling, hardcoded

m = 1; L = 10; %56.5962*scale;
gravity = 9.81;
[f, g, h] = getSystem11(nFterms, m, L);
fprintf('Running Example 11\n')

fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

%  Compute the polynomial approximations to the past future energy function
% [v] = approxPastEnergy(f, N, g, h, eta, degree, true);
[w] = pqr(f, g, h2q(h), eta, degree, true);

nX = 301; nY = nX;
xLim = pi; yLim = 5;
xPlot = linspace(-xLim, xLim, nX);
yPlot = linspace(-yLim, yLim, nY);
[X, Y] = meshgrid(xPlot, yPlot);

eFuture = zeros(nY, nX);
eFuture = zeros(nY, nX);

for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        eFuture(i, j) = 0.5 * kronPolyEval(w, x, degree);
        wRES(i, j) = computeResidualFutureHJB_2D_example11(gravity, L, g, h, eta, w, degree, x);
        if eFuture(i, j) < 0
            eFuture(i, j) = NaN;
        end
    end
end

fig1 = figure;
% ('Position',[600 50 1000 400])
% subplot(1,2,1)
contourf(X, Y, eFuture, 16, 'w'); hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
xticks([-pi, 0, pi])
xticklabels({'-\pi', '0', '\pi'})
%     axis equal
if degree > 2 && nFterms > 2
    %     caxis([0 1e4])
end

%         set(h, 'ylim', [0 1.5])
load(fullfile('utils', 'YlGnBuRescaled.mat'))
colormap(flip(YlGnBuRescaled))

if exportPlotData
    %     fprintf('Exporting matlab2tikz standalone tex file to: \n     plots/example11_futureEnergy_d%i_polynomial%i.tex\n',degree,nFterms)
    %     matlab2tikz('showInfo', false,'standalone',true,sprintf('plots/example11_futureEnergy_d%i_polynomial%i.tex',degree,nFterms))
    %     data = [ X(:) Y(:) eFuture(:) ];
    %     save plots/P.dat data -ASCII
    caxis([0 4e4])
    axis off
    fprintf('Exporting figure to: \n     plots/example11_futureEnergy_d%i_polynomial%i.pdf\n', degree, nFterms)
    exportgraphics(fig1, sprintf('plots/example11_futureEnergy_d%i_polynomial%i.pdf', degree, nFterms), 'ContentType', 'vector', 'BackgroundColor', 'none');
end
colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex', 'XTick', 0:10000:4e4, 'XTickLabel', {'0', '1e4', '2e4', '3e4', '4e4'});

% title('Future Energy Function')

fig2 = figure;
%     subplot(1,2,2)
pcolor(X, Y, log10(abs(wRES))); shading interp;
%     pcolor(X, Y, abs(wRES)); shading interp; colorbar;
% contourf(X, Y, abs(wRES), 50, 'w'); hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
xticks([-pi, 0, pi])
xticklabels({'-\pi', '0', '\pi'})
load('utils\YlGnBuRescaled.mat')
colormap(flip(YlGnBuRescaled))
% caxis([0 1e3])
caxis([-3 9])
% XTickLabels = cellstr(num2str(round(log10(XTick(:))), '10^%d'));

if exportPlotData
    %     fprintf('Exporting matlab2tikz standalone tex file to: \n     plots/example11_futureEnergy-HJB-Error_d%i_polynomial%i.tex\n',degree,nFterms)
    %     matlab2tikz('showInfo', false,'standalone',true,sprintf('plots/example11_futureEnergy-HJB-Error_d%i_polynomial%i.tex',degree,nFterms))
    axis off

    fprintf('Exporting figure to: \n     plots/example11_futureEnergy-HJB-Error_d%i_polynomial%i.pdf\n', degree, nFterms)
    exportgraphics(fig2, sprintf('plots/example11_futureEnergy-HJB-Error_d%i_polynomial%i.pdf', degree, nFterms), 'ContentType', 'vector', 'BackgroundColor', 'none');
end
colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex', 'XTick', -3:3:9, 'XTickLabel', {'1e-3', '1e0', '1e3', '1e6', '1e9'});
% title('HJB Residual')

end

function [res] = computeResidualFutureHJB_2D_example11(gravity, L, g, h, eta, w, degree, x)

w = w(1:degree);

%         constant B input
res = (0.5 * kronPolyDerivEval(w, x)) * [x(2); 3 * gravity / (2 * L) * sin(x(1))] ...
    - eta / 2 * 0.25 * kronPolyDerivEval(w, x) * g{1} * g{1}.' * kronPolyDerivEval(w, x).' ...
    + 0.5 * kronPolyEval(h, x).' * kronPolyEval(h, x);

%         Polynomial input
% res = (0.5 * kronPolyDerivEval(w, x)) * kronPolyEval(f, x) ...
%     - eta / 2 * 0.25 * kronPolyDerivEval(w, x) * (g{1} + kronPolyEval(g(2:end), x)) * (g{1} + kronPolyEval(g(2:end), x)).' * kronPolyDerivEval(w, x).' ...
%     + 0.5 * kronPolyEval(h, x).' * kronPolyEval(h, x);

end
