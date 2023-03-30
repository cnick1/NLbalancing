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
%              [2] Y. Kawano and J. M. A. Scherpen, “Model reduction by
%        differential balancing based on nonlinear hankel operators,”
%        IEEE Transactions on Automatic Control, vol. 62, no. 7,
%        pp. 3293–3308, Jul. 2017, doi: 10.1109/tac.2016.2628201.
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
[A, ~, C, N, g, ~, ~] = getSystem2(kawanoModel);
g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.
% n = size(A, 1);

eta = 0; % values should be between -\infty and 1.
% eta = 0.1; % values should be between -\infty and 1.
% eta=0.1 corresponds to gamma= 1.0541...
% since eta = 1 - 1/gamma^2;

% approximate the energy functions
[v] = approxPastEnergy(A, N, g(1:numGTermsApprox), C, eta, degree);
[w] = approxFutureEnergy(A, N, g(1:numGTermsApprox), C, eta, degree);

% compute the input-normal transformation approximation
[sigma, T] = inputNormalTransformation(v, w, degree - 1, true);

%% Plot the past and future energy functions
if (plotEnergy || plotBalancing)
    %  Plot the past and future energy functions in a neighborhood of the origin,
    %  first in the original coordinates, then in input normal
    nX = 101; nY = 101;
    xPlot = linspace(zmin, zmax, nX);
    yPlot = linspace(zmin, zmax, nY);
    [X, Y] = meshgrid(xPlot, yPlot);

    ePastOriginal = zeros(nX, nY);
    ePastInputNormalIdeal = zeros(nX, nY);
    ePastInputNormal = zeros(nX, nY);

    eFutureOriginal = zeros(nX, nY);
    eFutureInputNormal = zeros(nX, nY);

    for i = 1:nY
        for j = 1:nX
            z = [X(i, j); Y(i, j)];
            % First in original coordinates...
            ePastOriginal(i, j) = 0.5 * kronPolyEval(v, z, degree);
            eFutureOriginal(i, j) = 0.5 * kronPolyEval(w, z, degree);
            % ...then in input-normal coordinates.
            x = kronPolyEval(T, z, balancingDegree);
            ePastInputNormal(i, j) = 0.5 * kronPolyEval(v, x, degree);
            eFutureInputNormal(i, j) = 0.5 * kronPolyEval(w, x, degree);
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
        surf(X, Y, ePastOriginal)
        xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^-_\gamma(${\boldmath$x$}$)$', 'interpreter', 'latex');
        f1 = figure('Name', sprintf('Past energy function in original coordinates (degree %d approximation)', degree))
        contourf(X, Y, ePastOriginal); hold on;
        logMaxEPast = log10(max(max(ePastOriginal)));
        contour(X, Y, ePastOriginal, [0, logspace(-2, ceil(logMaxEPast), 20)] ./ (10 ^ (ceil(logMaxEPast) - logMaxEPast)))
        xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^-_\gamma(${\boldmath$x$}$)$', 'interpreter', 'latex');
        axis equal; colorbar

        figure('Name', sprintf('Past energy function in input-normal coordinates using degree %d transformation', balancingDegree))
        surf(X, Y, ePastInputNormal)
        xlabel('$z_1$', 'interpreter', 'latex'); ylabel('$z_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^-_\gamma(\Phi(${\boldmath$z$}$))$', 'interpreter', 'latex');
        f2 = figure('Name', sprintf('Past energy function in input-normal coordinates using degree %d transformation', balancingDegree))
        contourf(X, Y, ePastInputNormal); hold on;
        logMaxEPast = log10(max(max(ePastInputNormal)));
        contour(X, Y, ePastInputNormal, [0, logspace(-2, ceil(logMaxEPast), 20)] ./ (10 ^ (ceil(logMaxEPast) - logMaxEPast)))
        contour(X, Y, 0.5 * (X .^ 2 + Y .^ 2), 'k:', 'LineWidth', 2)
        xlabel('$z_1$', 'interpreter', 'latex'); ylabel('$z_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^-_\gamma(\Phi(${\boldmath$z$}$))$', 'interpreter', 'latex');
        axis equal; colorbar

        figure('Name', sprintf('Future energy function in original coordinates (degree %d approximation)', degree))
        surf(X, Y, eFutureOriginal)
        xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^+_\gamma(${\boldmath$x$}$)$', 'interpreter', 'latex');
        f3 = figure('Name', sprintf('Future energy function in original coordinates (degree %d approximation)', degree))
        contourf(X, Y, eFutureOriginal); hold on;
        logMaxEFuture = log10(max(max(eFutureOriginal)));
        contour(X, Y, eFutureOriginal, [0, logspace(-3, ceil(logMaxEFuture), 20)] ./ (10 ^ (ceil(logMaxEFuture) - logMaxEFuture)))
        xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^+_\gamma(${\boldmath$x$}$)$', 'interpreter', 'latex');
        axis equal; colorbar

        figure('Name', sprintf('Future energy function in input-normal coordinates using degree %d transformation', balancingDegree))
        surf(X, Y, eFutureInputNormal)
        xlabel('$z_1$', 'interpreter', 'latex'); ylabel('$z_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^+_\gamma(\Phi(${\boldmath$z$}$))$', 'interpreter', 'latex');
        f4 = figure('Name', sprintf('Future energy function in input-normal coordinates using degree %d transformation', balancingDegree))
        contourf(X, Y, eFutureInputNormal); hold on;
        logMaxEFuture = log10(max(max(eFutureInputNormal)))
        contour(X, Y, eFutureInputNormal, [0, logspace(-3, ceil(logMaxEFuture), 20)] ./ (10 ^ (ceil(logMaxEFuture) - logMaxEFuture)))
        xlabel('$z_1$', 'interpreter', 'latex'); ylabel('$z_2$', 'interpreter', 'latex');
        zlabel('$\mathcal{E}^+_\gamma(\Phi(${\boldmath$z$}$))$', 'interpreter', 'latex');
        axis equal; colorbar
    end
end

%% Approximate the singular value functions using Algorithm 2.
[c] = approximateSingularValueFunctions(T, w, sigma, degree - 2);

%% Generate data for plots of singular value functions
zRange = linspace(-0.2, 0.2, 51);
[svPlot, svSurface, xi] = plotSingularValueFunctions(sigma, c, zRange);
plotDifferentialSingularValueFunctions();

% TODO: Add export plot
%% Export plots
if exportPlotData
    figure(f1); exportgraphics(gca, 'plots/example2_pastEnergy_original.pdf', 'ContentType', 'vector');
    figure(f2); exportgraphics(gca, sprintf('plots/example2_pastEnergy_inputNormal_deg%d.pdf', balancingDegree), 'ContentType', 'vector');
    figure(f3); exportgraphics(gca, 'plots/example2_futureEnergy_original.pdf', 'ContentType', 'vector');
    figure(f4); exportgraphics(gca, sprintf('plots/example2_futureEnergy_inputNormal_deg%d.pdf', balancingDegree), 'ContentType', 'vector');
    
    figure(svPlot); exportgraphics(gca, 'plots/example2_singularValueFunctions.pdf', 'ContentType', 'vector');
    figure(svSurface); exportgraphics(gca, 'plots/example2_singularValueFunctionsSurface.pdf', 'ContentType', 'vector');
    
    fid = fopen('plots/ex2_singularValueFunctions.txt', 'w');
    fprintf(fid, '%g %g %g\n', [zRange; xi]);
    fclose(fid);
end
end
