function [v, w] = runExample2(degree, plotEnergy, plotBalancing, balancingDegree, numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel)
%runExample2 Runs the 2D example to plot energy functions as contour plots
%
%   Usage:  [v,w] = runExample2(degree,plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel)
%
%   runExample2() runs the default case of a quadratic model from [1] which
%                 is based on a model from [2]. The kawanoModel version is used in [3].
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
%   Reference: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki, “Nonlinear
%               balanced truncation: Part 1—computing energy functions,” arXiv,
%               Dec. 2022. doi: 10.48550/ARXIV.2209.07645
%              [2] Y. Kawano and J. M. A. Scherpen, “Model reduction by
%               differential balancing based on nonlinear hankel operators,”
%               IEEE Transactions on Automatic Control, vol. 62, no. 7,
%               pp. 3293–3308, Jul. 2017, doi: 10.1109/tac.2016.2628201.
%              [3] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” 2023.
%%

%% Process inputs
if nargin < 8
    if nargin < 7
        if nargin < 6
            if nargin < 5
                if nargin < 4
                    if nargin < 3
                        if nargin < 2
                            if nargin < 1
                                degree = 6;
                            end
                            plotEnergy = true;
                        end
                        plotBalancing = false;
                    end
                    balancingDegree = 3;
                end
                numGTermsModel = 1;
            end
            numGTermsApprox = numGTermsModel;
        end
        exportPlotData = false;
    end
    kawanoModel = false;
end

if plotBalancing
    dataRange = 0.2; %6.0
else
    dataRange = 1; %0.75;
end

%% Get model and compute energy functions
[f, g, h] = getSystem2(kawanoModel);
g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.
fprintf('Running Example 2\n')

if kawanoModel
    eta = 0; % values should be between -\infty and 1.
else
    eta = 0.1; % values should be between -\infty and 1.
end
% eta=0.1 corresponds to gamma= 1.0541...
% since eta = 1 - 1/gamma^2;

fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

%  Compute the polynomial approximations to the past future energy function
[v] = approxPastEnergy(f, g(1:numGTermsApprox), h, eta, degree, true);
[w] = approxFutureEnergy(f, g(1:numGTermsApprox), h, eta, degree, true);

%% Plot the past and future energy functions
if (plotEnergy || plotBalancing)
    nX = 301; nY = nX;
    xPlot = linspace(-dataRange, dataRange, nX);
    yPlot = linspace(-dataRange, dataRange, nY);
    ePast = zeros(nY, nX);
    eFuture = zeros(nY, nX);
    [X, Y] = meshgrid(xPlot, yPlot);

    for i = 1:nY
        for j = 1:nX
            x = [X(i, j); Y(i, j)];
            ePast(i, j) = 0.5 * kronPolyEval(v, x, degree);
            eFuture(i, j) = 0.5 * kronPolyEval(w, x, degree);
        end
    end
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaulttextinterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    fig1 = figure;
    contourf(X, Y, ePast, 16, 'w'); hold on;
    logMaxEPast = log10(max(max(ePast)));
    contour(X, Y, ePast, [0, logspace(-2, ceil(logMaxEPast), 20)] ./ (10 ^ (ceil(logMaxEPast) - logMaxEPast)))
    %    mesh(X,Y,ePast)
    xlabel('$x_1$', 'interpreter', 'latex');
    ylabel('$x_2$', 'interpreter', 'latex');
    h = colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
    set(gca, 'FontSize', 16)
    xticks([-1:1])
    yticks([-1:1])
    axis equal
    if kawanoModel
        caxis([0 80])
        set(h, 'ylim', [0 80])
        load(fullfile('utils', 'YlGnBuRescaled.mat'))
        colormap(flip(YlGnBuRescaled))
    end

    fig2 = figure
    contourf(X, Y, eFuture, 16, 'w'); hold on;
    logMaxEFuture = log10(max(max(eFuture)));
    contour(X, Y, eFuture, [0, logspace(-3, ceil(logMaxEFuture), 20)] ./ (10 ^ (ceil(logMaxEFuture) - logMaxEFuture)))
    xlabel('$x_1$', 'interpreter', 'latex');
    ylabel('$x_2$', 'interpreter', 'latex');
    h = colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
    set(gca, 'FontSize', 16)
    xticks([-1:1])
    yticks([-1:1])
    axis equal
    if kawanoModel
        caxis([0 1.5])
        set(h, 'ylim', [0 1.5])
        load(fullfile('utils', 'YlGnBuRescaled.mat'))
        colormap(flip(YlGnBuRescaled))
    end

    % Draw square around middle .2
    %     x1 =- .2;
    %     x2 = .2;
    %     y1 =- .2;
    %     y2 = .2;
    %     x = [x1, x2, x2, x1, x1];
    %     y = [y1, y1, y2, y2, y1];
    %     plot(x, y, 'w-');

    if exportPlotData
        % save('Ex2_RawData.mat', 'v', 'w')

        fid = fopen('plots/ex2_past_future.txt', 'w');
        fprintf(fid, '%g %g %g %g\n', [X(:), Y(:), ePast(:), eFuture(:)]);
        fclose(fid);

        exportgraphics(fig1, 'plots/PEF_p0_1.pdf', 'ContentType', 'vector');
        exportgraphics(fig2, 'plots/FEF_p0_1.pdf', 'ContentType', 'vector');
    end
    figure(fig1); title('Past Energy Function')
    figure(fig2); title('Future Energy Function')
end

%% Plot something about balancing(?)
if (plotBalancing)
    [sigma, T] = inputNormalTransformation(v, w, balancingDegree);
    nPts = 201;
    s = linspace(-2, 2, nPts);
    lin = T{1}(:, 1) * s;

    coord = lin;
    for k = 2:balancingDegree
        coord = coord + T{k}(:, 1) * s .^ k;
    end

    idxLin = zeros(1, nPts);
    linCount = 0;
    for i = 1:nPts
        if (norm(lin(:, i), inf) < dataRange)
            linCount = linCount + 1;
            idxLin(linCount) = i;
        end
    end
    idxLin = idxLin(1:linCount);

    idxCoord = zeros(1, nPts);
    coordCount = 0;
    for i = 1:nPts
        if (norm(coord(:, i), inf) < dataRange)
            coordCount = coordCount + 1;
            idxCoord(coordCount) = i;
        end
    end
    idxCoord = idxCoord(1:coordCount);

    figure(1); hold on
    plot(lin(1, idxLin), lin(2, idxLin), 'w+')
    plot(coord(1, idxCoord), coord(2, idxCoord), 'r+')

    figure(2); hold on
    plot(lin(1, idxLin), lin(2, idxLin), 'w+')
    plot(coord(1, idxCoord), coord(2, idxCoord), 'r+')
end
end
