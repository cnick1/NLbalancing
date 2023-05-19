function runExample2_regionOfAccuracy_res(exportData, varargin)
%runExample1_regionOfAccuracy Runs 1D ODE example to compare computed and
%analytical energy functions. This function plots a) error vs region size
%comparing degree of approximation for a polynomial function, and b) error
%vs region size comparing degree of assuemed model.
%
%   Usage:  runExample1_regionOfAccuracy()
%
%   Part of the NLbalancing repository.
%%

if nargin < 1
    exportData = false; %change
end

%% 1st Figure: all energy functions, big mess but just for me.

eta = 0; % values should be between -\infty and 1.

[A, ~, C, N, f, g, h] = getSystem2(true);

%  Compute the polynomial approximations to the future energy function
N=301;
xPlot = linspace(-1, 1, N);
yPlot = linspace(-1, 1, N);
[X, Y] = meshgrid(xPlot, yPlot);
[w] = approxFutureEnergy(f, N, g, h, eta, 8);

for d=2:8
    
    RES = computeResidualFutureHJB(f, g, h, eta, w, d);
    
    figure
    contourf(X,Y,abs(RES),16,'w'); colorbar;
end

end
