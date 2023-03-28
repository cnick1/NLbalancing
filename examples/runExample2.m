function [v, w] = runExample2(degree, plotEnergy, plotBalancing, balancingDegree, numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel, varargin)
%runExample2 Runs the 2D example to plot energy functions as contour plots
%
%   Usage:  [v,w] = runExample2(degree,plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel)
%
%   runExample2() runs the default case of a quadratic model from [1] which
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
    kawanoModel = 0;
end

if plotBalancing
    dataRange = 0.2; %6.0
else
    dataRange = 1; %0.75;
end

%% Get model and compute energy functions
[A, ~, C, N, g, ~, ~] = getSystem2(kawanoModel);
g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.
fprintf('Running Example 2\n')

eta = 0.1; % values should be between -\infty and 1.
% eta=0.1 corresponds to gamma= 1.0541...
% since eta = 1 - 1/gamma^2;

fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

%  Compute the polynomial approximations to the past future energy function
[v] = approxPastEnergy(A, N, g(1:numGTermsApprox), C, eta, degree, true);
[w] = approxFutureEnergy(A, N, g(1:numGTermsApprox), C, eta, degree, true);

wT{degree} = [];
vT{degree} = [];
for k = 2:degree
    wT{k} = w{k}.';
    vT{k} = v{k}.';
end

%% Plot the past and future energy functions
if (plotEnergy || plotBalancing)
    nX = 101; nY = 102;
    xPlot = linspace(-dataRange, dataRange, nX);
    yPlot = linspace(-dataRange, dataRange, nY);
    ePast = zeros(nY, nX);
    eFuture = zeros(nY, nX);
    [X, Y] = meshgrid(xPlot, yPlot);

    for i = 1:nY
        for j = 1:nX
            x = [X(i, j); Y(i, j)];
            %         WBar = vbar(wT,x);
            %         eFuture(i,j) = x.'*WBar*x;
            ePast(i, j) = 0.5 * kronPolyEval(vT, x, degree);
            eFuture(i, j) = 0.5 * kronPolyEval(wT, x, degree);
        end
    end
    figure
    contourf(X, Y, ePast)
    %    mesh(X,Y,ePast)
    xlabel('$x_1$', 'interpreter', 'latex');
    ylabel('$x_2$', 'interpreter', 'latex');
    colorbar('FontSize', 16)
    set(gca, 'FontSize', 20)

    figure
    contourf(X, Y, eFuture)
    xlabel('$x_1$', 'interpreter', 'latex');
    ylabel('$x_2$', 'interpreter', 'latex');
    colorbar('FontSize', 16)
    set(gca, 'FontSize', 20)

    if exportPlotData
        fid = fopen('plots/ex2_past_future.txt', 'w');
        fprintf(fid, '%g %g %g %g\n', [X(:); Y(:); ePast(:); eFuture(:)]);
        fclose(fid);
        % save('Ex2_RawData.mat', 'v', 'w')
        figure(1)
        exportgraphics(gca, 'plots/PEF_p0_1.pdf', 'ContentType', 'vector');

        figure(2)
        exportgraphics(gca, 'plots/FEF_p0_1.pdf', 'ContentType', 'vector');
    end
    figure(1); title('Past Energy Function')
    figure(2); title('Future Energy Function')
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
