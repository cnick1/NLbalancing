function runExample1_regionOfAccuracy_res(exportData)
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

xd = linspace(-6, 6, 250);
eta = 0.5; % values should be between -\infty and 1.
% eta=0.5 corresponds to gamma= sqrt(2)
% since eta = 1 - 1/gamma^2;

[f, g, h] = getSystem1();

%  Compute the polynomial approximations to the future energy function
d = 8;
[v] = approxPastEnergy(f, g, h, eta, d);

RES2 = computeResidualPastHJB(f, g, h, eta, v, 2, 6, 250);
RES4 = computeResidualPastHJB(f, g, h, eta, v, 4, 6, 250);
RES6 = computeResidualPastHJB(f, g, h, eta, v, 6, 6, 250);
RES8 = computeResidualPastHJB(f, g, h, eta, v, 8, 6, 250);

figure(3); colororder({'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30'});
set(gca, 'DefaultLineLineWidth', 2);
set(gca, 'DefaultAxesLineStyleOrder', {'-', '--', ':'});

% g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.

figure(3); plot(xd(1:10:end), 0 * xd(1:10:end), '+', 'LineWidth', 1); hold on;
plot(xd, abs(RES2), ...
    xd, abs(RES4), ...
    xd, abs(RES6), ...
    xd, abs(RES8), ...
    'LineWidth', 2)

legend('analytic', 'degree 2', 'degree 4', 'degree 6', 'degree 8', 'Location', 'southeast')
xlabel('$x$', 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold')
ylabel('$\mathcal{E}_\gamma^+$', 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold')
xlim([-6, 6]); ylim([-0.1, 3]);

intervalSize = xd(126:end);
Ep2_errorFun = zeros(length(intervalSize), 1);
Ep4_errorFun = zeros(length(intervalSize), 1);
Ep6_errorFun = zeros(length(intervalSize), 1);
Ep8_errorFun = zeros(length(intervalSize), 1);
for idx = 0:124
    Ep2_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), RES2(125 - idx:126 + idx) .^ 2)) / (2 * xd(126 + idx));
    Ep4_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), RES4(125 - idx:126 + idx) .^ 2)) / (2 * xd(126 + idx));
    Ep6_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), RES6(125 - idx:126 + idx) .^ 2)) / (2 * xd(126 + idx));
    Ep8_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), RES8(125 - idx:126 + idx) .^ 2)) / (2 * xd(126 + idx));
end

figure(4); colororder({'#D95319', '#EDB120', '#7E2F8E', '#77AC30'});
set(figure(4), 'DefaultLineLineWidth', 2);
set(figure(4), 'DefaultAxesLineStyleOrder', {'-', '--', ':'});

semilogy(intervalSize, Ep2_errorFun)
hold on;
semilogy(intervalSize, Ep4_errorFun)
semilogy(intervalSize, Ep6_errorFun)
semilogy(intervalSize, Ep8_errorFun)

if exportData
    fileName = sprintf('plots/example1_regionOfAccuracy_residual.dat');
    fprintf("Writing data to " + fileName + '\n')
    fileID = fopen(fileName, 'w');

    %print the header
    fprintf(fileID, 'intervalSize      & Ep2_errorFun    & Ep4_errorFun    & Ep6_errorFun    & Ep8_errorFun        \n');
    for i = 1:length(intervalSize)
        fprintf(fileID, '%12.6e      & ', intervalSize(i));
        fprintf(fileID, '%12.6e    & %12.6e    & %12.6e    & %12.6e       \n', Ep2_errorFun(i), Ep4_errorFun(i), Ep6_errorFun(i), Ep8_errorFun(i));
    end
    fclose(fileID);
    type(fileName);
end

end
