function [v, w] = runExample2_singularValueFunctions(degree, plotEnergy, plotBalancing, balancingDegree, numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel)
%runExample2_singularValueFunctions Runs the 2D example to plot energy
%functions as surfaces in the original and input normal coordinates. This
%function also plots the singular value functions.
%
%   Usage:  [v,w] = runExample2_singularValueFunctions(degree,plotEnergy,plotBalancing,balancingDegree
%                                   , numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel)
%
%   runExample2_singularValueFunctions() runs the default case of a quadratic model from [1] which
%                 is based on a model from [2].
%
%   Inputs:
%       degree          is the degree of energy function approximations
%       plotEnergy      is a logical variable to determine if a plot is made.
%
%       plotBalancing   is a logical variable to determine if a plot is made.
%       balancingDegree is a small integer setting the degree of approximation
%                         in the balancing transformation. Must be < degree.
%                         (used if plotBalancing = true)
%       numGTermsModel   Number of terms in the full order model
%       numGTermsApprox  Number of terms assumed when computing energy functions
%       exportPlotData   Boolean variable to determine if plots/data are exported
%       kawanoModel      Boolean variable modifying the output equation for
%                          the model to match [2]
%
%   Outputs:
%       v,w              are coefficients of the past and future energy
%                        function approximations, respectively.
%
%   The value of eta is set below.
%
%   References: [1] Nonlinear Balanced Truncation Model Reduction:
%        Part 1-Computing Energy Functions, by Kramer, Gugercin, and Borggaard.
%        arXiv:2209.07645.
%
%              [2] Y. Kawano and J. M. A. Scherpen, “Model reduction by
%        differential balancing based on nonlinear hankel operators,”
%        IEEE Transactions on Automatic Control, vol. 62, no. 7,
%        pp. 3293–3308, Jul. 2017, doi: 10.1109/tac.2016.2628201.
%
%               [3]"Scalable Computation of ℋ∞ Energy Functions for
%        Polynomial Control-Affine Systems", N. Corbin and B. Kramer
%        arXiv:
%
%   Part of the NLbalancing repository.
%% Process inputs
if nargin < 8
    if nargin < 7
        if nargin < 6
            if nargin < 5
                if nargin < 4
                    if nargin < 3
                        if nargin < 2
                            if nargin < 1
                                degree = 8;
                            end
                            plotEnergy = true;
                        end
                        plotBalancing = true;
                    end
                    balancingDegree = 2;
                end
                numGTermsModel = 1;
            end
            numGTermsApprox = numGTermsModel;
        end
        exportPlotData = false;
    end
    kawanoModel = false;
end

validateInputNormalPastEnergy = true;

zmax = 0.2; zmin = -zmax;
% Npts = 51;

%% Get model and compute energy functions and input normal transformation
[f, g, h] = getSystem2(kawanoModel);
g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.
% n = size(A, 1);

eta = 0; % values should be between -\infty and 1.
% eta = 0.1; % values should be between -\infty and 1.
% eta=0.1 corresponds to gamma= 1.0541...
% since eta = 1 - 1/gamma^2;

% approximate the energy functions
[v] = approxPastEnergy(f, g(1:numGTermsApprox), C, eta, degree);
[w] = approxFutureEnergy(f, g(1:numGTermsApprox), C, eta, degree);

% compute the input-normal transformation approximation
[sigma, T] = inputNormalTransformation(v, w, degree - 1, true);

%% Plot the past and future energy functions
if (plotEnergy || plotBalancing)
    %  Plot the past and future energy functions in a neighborhood of the origin,
    %  first in the original coordinates, then in input normal
    nX = 101; nY = 101;
    xPlot = linspace(zmin, zmax, nX);
    yPlot = linspace(zmin, zmax, nY);
    [Z1, Z2] = meshgrid(xPlot, yPlot);
    [X1, X2] = meshgrid(xPlot, yPlot);

    ePastOriginal = zeros(nX, nY);
    ePastInputNormalIdeal = zeros(nX, nY);
    ePastInputNormal = zeros(nX, nY);

    eFutureOriginal = zeros(nX, nY);
    eFutureInputNormal = zeros(nX, nY);

    for i = 1:nY
        for j = 1:nX
            x = [X1(i, j); X2(i, j)];
            % First in original coordinates...
            ePastOriginal(i, j) = 0.5 * kronPolyEval(v, x, degree);
            eFutureOriginal(i, j) = 0.5 * kronPolyEval(w, x, degree);
            % ...then in input-normal coordinates.
            z = [Z1(i, j); Z2(i, j)];
            xTrans = kronPolyEval(T, z, balancingDegree);
            ePastInputNormal(i, j) = 0.5 * kronPolyEval(v, xTrans, degree);
            eFutureInputNormal(i, j) = 0.5 * kronPolyEval(w, xTrans, degree);
            % compute ideal for checking errors
            ePastInputNormalIdeal(i, j) = 0.5 * (z.' * z);
        end
    end

    if (validateInputNormalPastEnergy)
        EminusError = max(max(abs(ePastInputNormal - ePastInputNormalIdeal)));
        fprintf('Transformed past energy function error at degree %d is %g\n', ...
            balancingDegree, EminusError);
        %[g,i] = max(max(abs(Eplot-Epast)))
    end

    if (plotEnergy)
        figure('Name', sprintf('Past energy function in original coordinates (degree %d approximation)', degree))
        surf(X1, X2, ePastOriginal)
        xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^-_\gamma(${\boldmath$x$}$)$', 'interpreter', 'latex');
        f1 = figure('Name', sprintf('Past energy function in original coordinates (degree %d approximation)', degree))
        contourf(X1, X2, ePastOriginal, 10, 'w:', 'LineWidth', 3); hold on;
        load(fullfile('utils', 'YlGnBuRescaled.mat'))
        colormap(flip(YlGnBuRescaled))
        caxis([0 80])
        %         logMaxEPast = log10(max(max(ePastOriginal)));
        %         contour(X1, X2, ePastOriginal, [0, logspace(-2, ceil(logMaxEPast), 20)] ./ (10 ^ (ceil(logMaxEPast) - logMaxEPast)))
        %         xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
        %         zlabel('$\mathcal{E}^-_\gamma(${\boldmath$x$}$)$', 'interpreter', 'latex');
        axis equal;
        %     colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex')
        %     set(gca, 'FontSize', 16)
        xticks([])
        yticks([])

        xPlot2 = linspace(- .25, .25, nX);
        yPlot2 = linspace(- .25, .25, nY);
        [Z12, Z22] = meshgrid(xPlot2, yPlot2);
        [X12, X22] = meshgrid(xPlot2, yPlot2);

        ePastOriginal2 = zeros(nX, nY);
        ePastInputNormalIdeal2 = zeros(nX, nY);
        ePastInputNormal2 = zeros(nX, nY);

        for i = 1:nY
            for j = 1:nX
                x = [X12(i, j); X22(i, j)];
                % First in original coordinates...
                ePastOriginal2(i, j) = 0.5 * kronPolyEval(v, x, degree);
                % ...then in input-normal coordinates.
                z = [Z12(i, j); Z22(i, j)];
                xTrans = kronPolyEval(T, z, balancingDegree);
                ePastInputNormal2(i, j) = 0.5 * kronPolyEval(v, xTrans, degree);
                eFutureInputNormal2(i, j) = 0.5 * kronPolyEval(w, xTrans, degree);
                % compute ideal for checking errors
                ePastInputNormalIdeal2(i, j) = 0.5 * (z.' * z);
            end
        end

        f2 = figure('Name', sprintf('Past energy function in input-normal coordinates using degree %d transformation', balancingDegree))
        contourf(Z12, Z22, ePastInputNormal2, 'w'); hold on;
        %         logMaxEPast = log10(max(max(ePastInputNormal)));
        %         contour(Z1, Z2, ePastInputNormal, [0, logspace(-2, ceil(logMaxEPast), 20)] ./ (10 ^ (ceil(logMaxEPast) - logMaxEPast)))
        contour(Z12, Z22, 0.5 * (Z12 .^ 2 + Z22 .^ 2), 'r:', 'LineWidth', 1.5)
        xlabel('$z_1$', 'interpreter', 'latex'); ylabel('$z_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^-_\gamma(\Phi(${\boldmath$z$}$))$', 'interpreter', 'latex');
        axis equal;
        xticks([- .2 0 .2])
        yticks([- .2 0 .2])
        if balancingDegree ~= 1
            colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex')
            ylabel('');
            yticks([]);
        end
        set(gca, 'FontSize', 16)
        load(fullfile('utils', 'YlGnBu.mat'))
        colormap(flip(YlGnBu))
        caxis([0 0.05])
        set(gca, 'Units', 'centimeters', 'Position', [3 3 6 6])
        %     , 'Units','centimeters', 'Position', [3 3 4 4])

        figure('Name', sprintf('Future energy function in original coordinates (degree %d approximation)', degree))
        surf(X1, X2, eFutureOriginal)
        xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^+_\gamma(${\boldmath$x$}$)$', 'interpreter', 'latex');
        f3 = figure('Name', sprintf('Future energy function in original coordinates (degree %d approximation)', degree))
        contourf(X1, X2, eFutureOriginal, 10, 'w:', 'LineWidth', 3); hold on;
        load(fullfile('utils', 'YlGnBuRescaled.mat'))
        colormap(flip(YlGnBuRescaled))
        caxis([0 1.5])
        %         logMaxEFuture = log10(max(max(eFutureOriginal)));
        %         contour(X1, X2, eFutureOriginal, [0, logspace(-3, ceil(logMaxEFuture), 20)] ./ (10 ^ (ceil(logMaxEFuture) - logMaxEFuture)))
        %         xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
        %         zlabel('$\mathcal{E}^+_\gamma(${\boldmath$x$}$)$', 'interpreter', 'latex');
        axis equal;
        %     colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex')
        %     set(gca, 'FontSize', 16)
        xticks([])
        yticks([])

        figure('Name', sprintf('Future energy function in input-normal coordinates using degree %d transformation', balancingDegree))
        surf(Z1, Z2, eFutureInputNormal)
        xlabel('$z_1$', 'interpreter', 'latex'); ylabel('$z_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^+_\gamma(\Phi(${\boldmath$z$}$))$', 'interpreter', 'latex');
        f4 = figure('Name', sprintf('Future energy function in input-normal coordinates using degree %d transformation', balancingDegree))
        contourf(Z1, Z2, eFutureInputNormal); hold on;
        logMaxEFuture = log10(max(max(eFutureInputNormal)))
        %         contour(Z1, Z2, eFutureInputNormal, [0, logspace(-3, ceil(logMaxEFuture), 20)] ./ (10 ^ (ceil(logMaxEFuture) - logMaxEFuture)),'w')
        xlabel('$z_1$', 'interpreter', 'latex'); ylabel('$z_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^+_\gamma(\Phi(${\boldmath$z$}$))$', 'interpreter', 'latex');
        axis equal;
        colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex')
        set(gca, 'FontSize', 16)
        xticks([- .2 0 .2])
        yticks([- .2 0 .2])
    end
end

% xlabel('$x_1$', 'interpreter', 'latex');
% ylabel('$x_2$', 'interpreter', 'latex');
% colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex')
% set(gca, 'FontSize', 16)

%% Approximate the singular value functions using Algorithm 2.
[c] = approximateSingularValueFunctions(T, w, sigma, degree - 2);

%% Generate data for plots of singular value functions
zRange = linspace(-0.2, 0.2, 51);
[svPlot, svSurface, xi] = plotSingularValueFunctions(sigma, c, zRange);
if kawanoModel
    [ctrbEnergy, obsvEnergy] = plotDifferentialSingularValueFunctions(1);
    [ctrbEnergy_zoomed, obsvEnergy_zoomed] = plotDifferentialSingularValueFunctions(.2);
end

%% Export plots
if exportPlotData
    exportgraphics(f1, 'plots/example2_pastEnergy_original.pdf', 'ContentType', 'vector');
    exportgraphics(f2, sprintf('plots/example2_pastEnergy_inputNormal_deg%d.pdf', balancingDegree), 'ContentType', 'vector');
    exportgraphics(f3, 'plots/example2_futureEnergy_original.pdf', 'ContentType', 'vector');
    exportgraphics(f4, sprintf('plots/example2_futureEnergy_inputNormal_deg%d.pdf', balancingDegree), 'ContentType', 'vector');

    exportgraphics(svPlot, 'plots/example2_singularValueFunctions.pdf', 'ContentType', 'vector');
    exportgraphics(svSurface, 'plots/example2_singularValueFunctionsSurface.pdf', 'ContentType', 'vector');

    fid = fopen('plots/ex2_singularValueFunctions.txt', 'w');
    fprintf(fid, '%g %g %g\n', [zRange; xi]);
    fclose(fid);

    if kawanoModel
        exportgraphics(ctrbEnergy, 'plots/example2_kawanoCtrbEnergy.pdf', 'ContentType', 'vector');
        exportgraphics(obsvEnergy, 'plots/example2_kawanoObsvEnergy.pdf', 'ContentType', 'vector');
        exportgraphics(ctrbEnergy_zoomed, 'plots/example2_kawanoCtrbEnergy_zoomed.pdf', 'ContentType', 'vector');
        exportgraphics(obsvEnergy_zoomed, 'plots/example2_kawanoObsvEnergy_zoomed.pdf', 'ContentType', 'vector');
    end

end
end
