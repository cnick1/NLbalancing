function runExample1_regionOfAccuracy(exportData, varargin)
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

for numGTermsModel = 3:-1:1
    g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.

    %  Compute the polynomial approximations to the future energy function
    d = 8;
    [w] = approxFutureEnergy(A, N, g, C, eta, d);
    w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};
    Ef2 = 0.5 * w2 * xd .^ 2; Ef3 = Ef2 + 0.5 * w3 * xd .^ 3; Ef4 = Ef3 + 0.5 * w4 * xd .^ 4; Ef5 = Ef4 + 0.5 * w5 * xd .^ 5; Ef6 = Ef5 + 0.5 * w6 * xd .^ 6; Ef7 = Ef6 + 0.5 * w7 * xd .^ 7; Ef8 = Ef7 + 0.5 * w8 * xd .^ 8;

    figure(1); plot(xd(1:10:end), EPlusAnalytic(1:10:end), '+', 'LineWidth', 1); hold on;
    plot(xd, Ef2, ...
        xd, Ef4, ...
        xd, Ef6, ...
        xd, Ef8, ...
        'LineWidth', 2)
end

legend('analytic', 'degree 2', 'degree 4', 'degree 6', 'degree 8', 'Location', 'southeast')
xlabel('$x$', 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold')
ylabel('$\mathcal{E}_\gamma^+$', 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold')
xlim([-6, 6]); ylim([0, 10]);

%% 2nd Figure: all errors, big mess but just for me.

xd = linspace(-6, 6, 250);
eta = 0.5; % values should be between -\infty and 1.

[A, B, C, N, f, g, h] = getSystem1();

figure(2); colororder({'#D95319', '#EDB120', '#7E2F8E', '#77AC30'});
set(gca, 'DefaultLineLineWidth', 2);
set(gca, 'DefaultAxesLineStyleOrder', {'-', '--', ':'});

for numGTermsModel = 3:-1:1
    g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.

    %  Compute the polynomial approximations to the future energy function
    d = 8;
    [w] = approxFutureEnergy(A, N, g, C, eta, d);
    w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};
    Ef2 = 0.5 * w2 * xd .^ 2; Ef3 = Ef2 + 0.5 * w3 * xd .^ 3; Ef4 = Ef3 + 0.5 * w4 * xd .^ 4; Ef5 = Ef4 + 0.5 * w5 * xd .^ 5; Ef6 = Ef5 + 0.5 * w6 * xd .^ 6; Ef7 = Ef6 + 0.5 * w7 * xd .^ 7; Ef8 = Ef7 + 0.5 * w8 * xd .^ 8;

    Ef2_error = Ef2.' - EPlusAnalytic; Ef2_error = Ef2_error .^ 2;
    Ef4_error = Ef4.' - EPlusAnalytic; Ef4_error = Ef4_error .^ 2;
    Ef6_error = Ef6.' - EPlusAnalytic; Ef6_error = Ef6_error .^ 2;
    Ef8_error = Ef8.' - EPlusAnalytic; Ef8_error = Ef8_error .^ 2;

    intervalSize = xd(126:end);
    Ef2_errorFun = zeros(length(intervalSize), 1);
    Ef4_errorFun = zeros(length(intervalSize), 1);
    Ef6_errorFun = zeros(length(intervalSize), 1);
    Ef8_errorFun = zeros(length(intervalSize), 1);
    for idx = 0:124
        Ef2_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef2_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef4_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef4_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef6_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef6_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef8_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef8_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
    end

    semilogy(intervalSize, Ef2_errorFun)
    hold on;
    semilogy(intervalSize, Ef4_errorFun)
    semilogy(intervalSize, Ef6_errorFun)
    semilogy(intervalSize, Ef8_errorFun)
end
legend('degree 2', 'degree 4', 'degree 6', 'degree 8', 'Location', 'southeast')
xlabel('Region radius from origin', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
ylabel('L2 error for energy function approximation', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
grid on;

%% 3rd Figure: only full model errors

xd = linspace(-6, 6, 250);
eta = 0.5; % values should be between -\infty and 1.

[A, B, C, N, f, g, h] = getSystem1();

figure(3); colororder({'#D95319', '#EDB120', '#7E2F8E', '#77AC30'});
set(gca, 'DefaultLineLineWidth', 2);
set(gca, 'DefaultAxesLineStyleOrder', {'-', '--', ':'});

for numGTermsModel = 3
    g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.

    %  Compute the polynomial approximations to the future energy function
    d = 8;
    [w] = approxFutureEnergy(A, N, g, C, eta, d);
    w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};
    Ef2 = 0.5 * w2 * xd .^ 2; Ef3 = Ef2 + 0.5 * w3 * xd .^ 3; Ef4 = Ef3 + 0.5 * w4 * xd .^ 4; Ef5 = Ef4 + 0.5 * w5 * xd .^ 5; Ef6 = Ef5 + 0.5 * w6 * xd .^ 6; Ef7 = Ef6 + 0.5 * w7 * xd .^ 7; Ef8 = Ef7 + 0.5 * w8 * xd .^ 8;

    Ef2_error = Ef2.' - EPlusAnalytic; Ef2_error = Ef2_error .^ 2;
    Ef4_error = Ef4.' - EPlusAnalytic; Ef4_error = Ef4_error .^ 2;
    Ef6_error = Ef6.' - EPlusAnalytic; Ef6_error = Ef6_error .^ 2;
    Ef8_error = Ef8.' - EPlusAnalytic; Ef8_error = Ef8_error .^ 2;

    intervalSize = xd(126:end);
    Ef2_errorFun = zeros(length(intervalSize), 1);
    Ef4_errorFun = zeros(length(intervalSize), 1);
    Ef6_errorFun = zeros(length(intervalSize), 1);
    Ef8_errorFun = zeros(length(intervalSize), 1);
    for idx = 0:124
        Ef2_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef2_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef4_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef4_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef6_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef6_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef8_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef8_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
    end

    semilogy(intervalSize, Ef2_errorFun)
    hold on;
    semilogy(intervalSize, Ef4_errorFun)
    semilogy(intervalSize, Ef6_errorFun)
    semilogy(intervalSize, Ef8_errorFun)
end
legend('degree 2', 'degree 4', 'degree 6', 'degree 8', 'Location', 'southeast')
xlabel('Region radius from origin', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
ylabel('L2 error for energy function approximation', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
grid on;

Ef2_errorFun_full = Ef2_errorFun;
Ef4_errorFun_full = Ef4_errorFun;
Ef6_errorFun_full = Ef6_errorFun;
Ef8_errorFun_full = Ef8_errorFun;

%% 4th Figure: only QB model errors

xd = linspace(-6, 6, 250);
eta = 0.5; % values should be between -\infty and 1.

[A, B, C, N, f, g, h] = getSystem1();

figure(4); colororder({'#D95319', '#EDB120', '#7E2F8E', '#77AC30'});
set(gca, 'DefaultLineLineWidth', 2);

for numGTermsModel = 2
    g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.

    %  Compute the polynomial approximations to the future energy function
    d = 8;
    [w] = approxFutureEnergy(A, N, g, C, eta, d);
    w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};
    Ef2 = 0.5 * w2 * xd .^ 2; Ef3 = Ef2 + 0.5 * w3 * xd .^ 3; Ef4 = Ef3 + 0.5 * w4 * xd .^ 4; Ef5 = Ef4 + 0.5 * w5 * xd .^ 5; Ef6 = Ef5 + 0.5 * w6 * xd .^ 6; Ef7 = Ef6 + 0.5 * w7 * xd .^ 7; Ef8 = Ef7 + 0.5 * w8 * xd .^ 8;

    Ef2_error = Ef2.' - EPlusAnalytic; Ef2_error = Ef2_error .^ 2;
    Ef4_error = Ef4.' - EPlusAnalytic; Ef4_error = Ef4_error .^ 2;
    Ef6_error = Ef6.' - EPlusAnalytic; Ef6_error = Ef6_error .^ 2;
    Ef8_error = Ef8.' - EPlusAnalytic; Ef8_error = Ef8_error .^ 2;

    intervalSize = xd(126:end);
    Ef2_errorFun = zeros(length(intervalSize), 1);
    Ef4_errorFun = zeros(length(intervalSize), 1);
    Ef6_errorFun = zeros(length(intervalSize), 1);
    Ef8_errorFun = zeros(length(intervalSize), 1);
    for idx = 0:124
        Ef2_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef2_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef4_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef4_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef6_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef6_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef8_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef8_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
    end

    semilogy(intervalSize, Ef2_errorFun)
    hold on;
    semilogy(intervalSize, Ef4_errorFun, '--')
    semilogy(intervalSize, Ef6_errorFun, '--')
    semilogy(intervalSize, Ef8_errorFun, '--')
end
legend('degree 2', 'degree 4', 'degree 6', 'degree 8', 'Location', 'southeast')
xlabel('Region radius from origin', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
ylabel('L2 error for energy function approximation', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
grid on;

Ef2_errorFun_QB = Ef2_errorFun;
Ef4_errorFun_QB = Ef4_errorFun;
Ef6_errorFun_QB = Ef6_errorFun;
Ef8_errorFun_QB = Ef8_errorFun;

%% 5th Figure: only Q model errors

xd = linspace(-6, 6, 250);
eta = 0.5; % values should be between -\infty and 1.

[A, B, C, N, f, g, h] = getSystem1();

figure(5); colororder({'#D95319', '#EDB120', '#7E2F8E', '#77AC30'});
set(gca, 'DefaultLineLineWidth', 2);

for numGTermsModel = 1
    g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.

    %  Compute the polynomial approximations to the future energy function
    d = 8;
    [w] = approxFutureEnergy(A, N, g, C, eta, d);
    w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};
    Ef2 = 0.5 * w2 * xd .^ 2; Ef3 = Ef2 + 0.5 * w3 * xd .^ 3; Ef4 = Ef3 + 0.5 * w4 * xd .^ 4; Ef5 = Ef4 + 0.5 * w5 * xd .^ 5; Ef6 = Ef5 + 0.5 * w6 * xd .^ 6; Ef7 = Ef6 + 0.5 * w7 * xd .^ 7; Ef8 = Ef7 + 0.5 * w8 * xd .^ 8;

    Ef2_error = Ef2.' - EPlusAnalytic; Ef2_error = Ef2_error .^ 2;
    Ef4_error = Ef4.' - EPlusAnalytic; Ef4_error = Ef4_error .^ 2;
    Ef6_error = Ef6.' - EPlusAnalytic; Ef6_error = Ef6_error .^ 2;
    Ef8_error = Ef8.' - EPlusAnalytic; Ef8_error = Ef8_error .^ 2;

    intervalSize = xd(126:end);
    Ef2_errorFun = zeros(length(intervalSize), 1);
    Ef4_errorFun = zeros(length(intervalSize), 1);
    Ef6_errorFun = zeros(length(intervalSize), 1);
    Ef8_errorFun = zeros(length(intervalSize), 1);
    for idx = 0:124
        Ef2_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef2_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef4_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef4_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef6_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef6_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
        Ef8_errorFun(idx + 1) = sqrt(trapz(xd(125 - idx:126 + idx), Ef8_error(125 - idx:126 + idx))) / (2 * xd(126 + idx));
    end

    semilogy(intervalSize, Ef2_errorFun)
    hold on;
    semilogy(intervalSize, Ef4_errorFun, ':')
    semilogy(intervalSize, Ef6_errorFun, ':')
    semilogy(intervalSize, Ef8_errorFun, ':')
end
legend('degree 2', 'degree 4', 'degree 6', 'degree 8', 'Location', 'southeast')
xlabel('Region radius from origin', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
ylabel('L2 error for energy function approximation', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
grid on;

Ef2_errorFun_Q = Ef2_errorFun;
Ef4_errorFun_Q = Ef4_errorFun;
Ef6_errorFun_Q = Ef6_errorFun;
Ef8_errorFun_Q = Ef8_errorFun;

if exportData
    fileName = sprintf('plots/example1_regionOfAccuracy.dat');
    fprintf("Writing data to " + fileName + '\n')
    fileID = fopen(fileName, 'w');

    %print the header
    fprintf(fileID, 'intervalSize      & Ef2_errorFun_full    & Ef4_errorFun_full    & Ef6_errorFun_full    & Ef8_errorFun_full    & Ef2_errorFun_QB    & Ef4_errorFun_QB    & Ef6_errorFun_QB    & Ef8_errorFun_QB    & Ef2_errorFun_Q    & Ef4_errorFun_Q    & Ef6_errorFun_Q    & Ef8_errorFun_Q        \n');
    for i = 1:length(intervalSize)
        fprintf(fileID, '%12.6e      & ', intervalSize(i));
        fprintf(fileID, '%12.6e         & %12.6e         & %12.6e         & %12.6e         & %12.6e       & %12.6e       & %12.6e       & %12.6e       & %12.6e      & %12.6e      & %12.6e      & %12.6e     \n', Ef2_errorFun_full(i), Ef4_errorFun_full(i), Ef6_errorFun_full(i), Ef8_errorFun_full(i), Ef2_errorFun_QB(i), Ef4_errorFun_QB(i), Ef6_errorFun_QB(i), Ef8_errorFun_QB(i), Ef2_errorFun_Q(i), Ef4_errorFun_Q(i), Ef6_errorFun_Q(i), Ef8_errorFun_Q(i));
    end
    fclose(fileID);
    type(fileName);
end


%% HJB Residual Error
xd = linspace(-6, 6, 301);
eta = 0.5; % values should be between -\infty and 1.
% eta=0.5 corresponds to gamma= sqrt(2)
% since eta = 1 - 1/gamma^2;

[A, B, C, N, f, g, h] = getSystem1();



end

function [E] = EgammaPlusOG(xd, A, B, C, N, g, eta)

syms x;

a = -eta / 2 * (g{1} + g{2} * x + g{3} * x ^ 2) ^ 2;
b = A * x + N * x ^ 2;
c = 1/2 * C ^ 2 * x ^ 2;

dEx1 = (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a);
dEx2 = (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a);

Ex1 = int(dEx1, x);
Ex2 = int(dEx2, x);

y = piecewise(x < 0, Ex2, x >= 0, Ex1);

fplot(y, [-6, 6], '-.+', 'LineWidth', 2)
hold on

E = y;
end

function [Ex] = EgammaPlusNumerical(xd, A, B, C, N, g, eta)

x = sym('x');

a = -eta / 2 * (g{1} + g{2} * x + g{3} * x ^ 2) ^ 2;
b = A * x + N * x ^ 2;
c = 1/2 * C ^ 2 * x ^ 2;

dEx2 = (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a);
dEx1 = (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a);

% Ex1 = int(dEx1, x);
% Ex2 = int(dEx2, x);

% y = piecewise(x < 0, Ex1, x >= 0, Ex2);

i = 1; Ex = zeros(length(xd), 1);
for xi = xd
    if xi < 0
        Ex(i) = vpaintegral(dEx1, x, 0, xi);
    else
        Ex(i) = vpaintegral(dEx2, x, 0, xi);
    end
    i = i + 1;
end

% fplot(y, [-6, 6], '-.+', 'LineWidth', 2)
% hold on
% plot(xd(1:10:end),Ex(1:10:end),'+', 'LineWidth', 1)

end
