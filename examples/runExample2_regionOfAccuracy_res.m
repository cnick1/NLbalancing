function runExample2_regionOfAccuracy_res()
%runExample2_regionOfAccuracy_res Runs 1D ODE example to compare computed and
%analytical energy functions. This function plots a) error vs region size
%comparing degree of approximation for a polynomial function, and b) error
%vs region size comparing degree of assuemed model.
%
%   Usage:  runExample2_regionOfAccuracy_res()
%
%   Part of the NLbalancing repository.
%%

%% 1st Figure: all energy functions, big mess but just for me.

eta = 0; % values should be between -\infty and 1.

[f, g, h] = getSystem2(true);

%  Compute the polynomial approximations to the future energy function
N = 301;
xPlot = linspace(-1, 1, N);
yPlot = linspace(-1, 1, N);
[X, Y] = meshgrid(xPlot, yPlot);
[v] = approxPastEnergy(f, g, h, eta, 4, true);
[w] = approxFutureEnergy(f, g, h, eta, 4, true);

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

for d = 4

    vRES = computeResidualPastHJB(f, g, h, eta, v, d, 1, 301);
    wRES = computeResidualFutureHJB(f, g, h, eta, w, d, 1, 301);

    fig1 = figure
    contourf(X, Y, abs(vRES), 16, 'w'); colorbar;
    xlabel('$x_1$', 'interpreter', 'latex');
    ylabel('$x_2$', 'interpreter', 'latex');
    h = colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
    set(gca, 'FontSize', 16)
    xticks([-1:1])
    yticks([-1:1])
    axis equal
    load('utils\YlGnBuRescaled.mat')
    colormap(flip(YlGnBuRescaled))
    fprintf('The residual of the HJB equation on the unit square is %g\n', norm(vRES, 'inf'));

    fig2 = figure
    pcolor(X, Y, abs(wRES)); shading interp; colorbar;
    xlabel('$x_1$', 'interpreter', 'latex');
    ylabel('$x_2$', 'interpreter', 'latex');
    h = colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
    set(gca, 'FontSize', 16)
    xticks([-1:1])
    yticks([-1:1])
    h = get(gca, 'DataAspectRatio')
    if h(3) == 1
        set(gca, 'DataAspectRatio', [1 1 1 / max(h(1:2))])
    else
        set(gca, 'DataAspectRatio', [1 1 h(3)])
    end
    load('utils\YlGnBuRescaled.mat')
    colormap(flip(YlGnBuRescaled))
    fprintf('The residual of the HJB equation on the unit square is %g\n', norm(wRES, 'inf'));

    exportgraphics(fig1, 'plots/example2_pastEnergy_residual.pdf', 'ContentType', 'vector');
    exportgraphics(fig2, 'plots/example2_futureEnergy_residual.pdf', 'ContentType', 'vector');

end

end
