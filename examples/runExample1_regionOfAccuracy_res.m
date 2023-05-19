function runExample1_regionOfAccuracy_res(exportData, varargin)
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

[A, B, C, N, f, g, h] = getSystem1();
EPlusAnalytic = EgammaPlusNumerical(xd, A, B, C, N, g, eta);

figure(1); colororder({'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30'});
set(gca, 'DefaultLineLineWidth', 2);
set(gca, 'DefaultAxesLineStyleOrder', {'-', '--', ':'});


% g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.

%  Compute the polynomial approximations to the future energy function
d = 8;
[w] = approxFutureEnergy(f, N, g, C, eta, d);
w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};
Ef2 = 0.5 * w2 * xd .^ 2; Ef3 = Ef2 + 0.5 * w3 * xd .^ 3; Ef4 = Ef3 + 0.5 * w4 * xd .^ 4; Ef5 = Ef4 + 0.5 * w5 * xd .^ 5; Ef6 = Ef5 + 0.5 * w6 * xd .^ 6; Ef7 = Ef6 + 0.5 * w7 * xd .^ 7; Ef8 = Ef7 + 0.5 * w8 * xd .^ 8;

figure(1); plot(xd(1:10:end), EPlusAnalytic(1:10:end), '+', 'LineWidth', 1); hold on;
plot(xd, Ef2, ...
    xd, Ef4, ...
    xd, Ef6, ...
    xd, Ef8, ...
    'LineWidth', 2)


legend('analytic', 'degree 2', 'degree 4', 'degree 6', 'degree 8', 'Location', 'southeast')
xlabel('$x$', 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold')
ylabel('$\mathcal{E}_\gamma^+$', 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold')
xlim([-6, 6]); ylim([0, 10]);


end

