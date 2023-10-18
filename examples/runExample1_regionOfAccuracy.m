function runExample1_regionOfAccuracy(exportData)
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
    exportData = false;
end

%% 3rd Figure: only full model errors

xd = linspace(-6, 6, 250);
eta = 0.5; % values should be between -\infty and 1.

[f, g, h] = getSystem1();
EMinusAnalytic = EgammaMinusNumerical(xd, f, g, h, eta);

%  Compute the polynomial approximations to the Past energy function
d = 8;
[v] = approxPastEnergy(f, g, h, eta, d);
v2 = v{2}; v3 = v{3}; v4 = v{4}; v5 = v{5}; v6 = v{6}; v7 = v{7}; v8 = v{8};
Ep2 = 0.5 * v2 * xd .^ 2; Ep3 = Ep2 + 0.5 * v3 * xd .^ 3; Ep4 = Ep3 + 0.5 * v4 * xd .^ 4; Ep5 = Ep4 + 0.5 * v5 * xd .^ 5; Ep6 = Ep5 + 0.5 * v6 * xd .^ 6; Ep7 = Ep6 + 0.5 * v7 * xd .^ 7; Ep8 = Ep7 + 0.5 * v8 * xd .^ 8;

% Linfty error

Ep2_error = abs(EMinusAnalytic - Ep2);
Ep4_error = abs(EMinusAnalytic - Ep4);
Ep6_error = abs(EMinusAnalytic - Ep6);
Ep8_error = abs(EMinusAnalytic - Ep8);

intervalSize = xd(126:end);
Ep2_errorFun = zeros(length(intervalSize), 1);
Ep4_errorFun = zeros(length(intervalSize), 1);
Ep6_errorFun = zeros(length(intervalSize), 1);
Ep8_errorFun = zeros(length(intervalSize), 1);
for idx = 0:124
    Ep2_errorFun(idx + 1) = max(Ep2_error(125 - idx:126 + idx));
    Ep4_errorFun(idx + 1) = max(Ep4_error(125 - idx:126 + idx));
    Ep6_errorFun(idx + 1) = max(Ep6_error(125 - idx:126 + idx));
    Ep8_errorFun(idx + 1) = max(Ep8_error(125 - idx:126 + idx));
end

Ep2_errorFun_full = Ep2_errorFun;
Ep4_errorFun_full = Ep4_errorFun;
Ep6_errorFun_full = Ep6_errorFun;
Ep8_errorFun_full = Ep8_errorFun;

%% 4th Figure: only QB model errors

xd = linspace(-6, 6, 250);
eta = 0.5; % values should be between -\infty and 1.

[f, g, h] = getSystem1();

numGTermsModel = 2;
g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.

%  Compute the polynomial approximations to the Past energy function
d = 8;
[v] = approxPastEnergy(f, g, h, eta, d);
v2 = v{2}; v3 = v{3}; v4 = v{4}; v5 = v{5}; v6 = v{6}; v7 = v{7}; v8 = v{8};
Ep2 = 0.5 * v2 * xd .^ 2; Ep3 = Ep2 + 0.5 * v3 * xd .^ 3; Ep4 = Ep3 + 0.5 * v4 * xd .^ 4; Ep5 = Ep4 + 0.5 * v5 * xd .^ 5; Ep6 = Ep5 + 0.5 * v6 * xd .^ 6; Ep7 = Ep6 + 0.5 * v7 * xd .^ 7; Ep8 = Ep7 + 0.5 * v8 * xd .^ 8;

% Linfty error

Ep2_error = abs(EMinusAnalytic - Ep2);
Ep4_error = abs(EMinusAnalytic - Ep4);
Ep6_error = abs(EMinusAnalytic - Ep6);
Ep8_error = abs(EMinusAnalytic - Ep8);

intervalSize = xd(126:end);
Ep2_errorFun = zeros(length(intervalSize), 1);
Ep4_errorFun = zeros(length(intervalSize), 1);
Ep6_errorFun = zeros(length(intervalSize), 1);
Ep8_errorFun = zeros(length(intervalSize), 1);
for idx = 0:124
    Ep2_errorFun(idx + 1) = max(Ep2_error(125 - idx:126 + idx));
    Ep4_errorFun(idx + 1) = max(Ep4_error(125 - idx:126 + idx));
    Ep6_errorFun(idx + 1) = max(Ep6_error(125 - idx:126 + idx));
    Ep8_errorFun(idx + 1) = max(Ep8_error(125 - idx:126 + idx));
end

Ep2_errorFun_QB = Ep2_errorFun;
Ep4_errorFun_QB = Ep4_errorFun;
Ep6_errorFun_QB = Ep6_errorFun;
Ep8_errorFun_QB = Ep8_errorFun;

%% 5th Figure: only Q model errors

xd = linspace(-6, 6, 250);
eta = 0.5; % values should be between -\infty and 1.

[f, g, h] = getSystem1();

numGTermsModel = 1;
g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.

%  Compute the polynomial approximations to the Past energy function
d = 8;
[v] = approxPastEnergy(f, g, h, eta, d);
v2 = v{2}; v3 = v{3}; v4 = v{4}; v5 = v{5}; v6 = v{6}; v7 = v{7}; v8 = v{8};
Ep2 = 0.5 * v2 * xd .^ 2; Ep3 = Ep2 + 0.5 * v3 * xd .^ 3; Ep4 = Ep3 + 0.5 * v4 * xd .^ 4; Ep5 = Ep4 + 0.5 * v5 * xd .^ 5; Ep6 = Ep5 + 0.5 * v6 * xd .^ 6; Ep7 = Ep6 + 0.5 * v7 * xd .^ 7; Ep8 = Ep7 + 0.5 * v8 * xd .^ 8;

% Linfty error

Ep2_error = abs(EMinusAnalytic - Ep2);
Ep4_error = abs(EMinusAnalytic - Ep4);
Ep6_error = abs(EMinusAnalytic - Ep6);
Ep8_error = abs(EMinusAnalytic - Ep8);

intervalSize = xd(126:end);
Ep2_errorFun = zeros(length(intervalSize), 1);
Ep4_errorFun = zeros(length(intervalSize), 1);
Ep6_errorFun = zeros(length(intervalSize), 1);
Ep8_errorFun = zeros(length(intervalSize), 1);
for idx = 0:124
    Ep2_errorFun(idx + 1) = max(Ep2_error(125 - idx:126 + idx));
    Ep4_errorFun(idx + 1) = max(Ep4_error(125 - idx:126 + idx));
    Ep6_errorFun(idx + 1) = max(Ep6_error(125 - idx:126 + idx));
    Ep8_errorFun(idx + 1) = max(Ep8_error(125 - idx:126 + idx));
end

Ep2_errorFun_Q = Ep2_errorFun;
Ep4_errorFun_Q = Ep4_errorFun;
Ep6_errorFun_Q = Ep6_errorFun;
Ep8_errorFun_Q = Ep8_errorFun;

% Plot the data for figure 3
figure(3); colororder({'#D95319', '#EDB120', '#7E2F8E', '#77AC30'});
set(figure(3), 'DefaultLineLineWidth', 2);
set(figure(3), 'DefaultAxesLineStyleOrder', {'-', '--', ':'});

semilogy(intervalSize, Ep2_errorFun_full)
hold on;
semilogy(intervalSize, Ep4_errorFun_full)
semilogy(intervalSize, Ep6_errorFun_full)
semilogy(intervalSize, Ep8_errorFun_full)

legend('degree 2', 'degree 4', 'degree 6', 'degree 8', 'Location', 'southeast')
xlabel('Region radius from origin', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
ylabel('Linfty error for energy function approximation', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
grid on;
ylim([1e-16 2e2])

% Plot the data for figure 4
figure(5); colororder({'#D95319', '#EDB120', '#7E2F8E', '#77AC30'});
set(figure(5), 'DefaultLineLineWidth', 2);

semilogy(intervalSize, Ep4_errorFun_full)
hold on;
semilogy(intervalSize, Ep8_errorFun_full, ':')
semilogy(intervalSize, Ep8_errorFun_QB, ':')
semilogy(intervalSize, Ep8_errorFun_Q, ':')

legend('degree 4', 'degree 8 full', 'degree 8 QB', 'degree 8 Q', 'Location', 'southeast')
xlabel('Region radius from origin', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
ylabel('Linfty error for energy function approximation', 'interpreter', 'latex', 'FontSize', 12, 'fontweight', 'bold')
grid on;
ylim([1e-16 2e2])

if exportData
    fileName = sprintf('plots/example1_regionOfAccuracy.dat');
    fprintf("Writing data to " + fileName + '\n')
    fileID = fopen(fileName, 'w');

    %print the header
    fprintf(fileID, 'intervalSize      & Ep2_errorFun_full    & Ep4_errorFun_full    & Ep6_errorFun_full    & Ep8_errorFun_full    & Ep2_errorFun_QB    & Ep4_errorFun_QB    & Ep6_errorFun_QB    & Ep8_errorFun_QB    & Ep2_errorFun_Q    & Ep4_errorFun_Q    & Ep6_errorFun_Q    & Ep8_errorFun_Q        \n');
    for i = 1:length(intervalSize)
        fprintf(fileID, '%12.6e      & ', intervalSize(i));
        fprintf(fileID, '%12.6e         & %12.6e         & %12.6e         & %12.6e         & %12.6e       & %12.6e       & %12.6e       & %12.6e       & %12.6e      & %12.6e      & %12.6e      & %12.6e     \n', Ep2_errorFun_full(i), Ep4_errorFun_full(i), Ep6_errorFun_full(i), Ep8_errorFun_full(i), Ep2_errorFun_QB(i), Ep4_errorFun_QB(i), Ep6_errorFun_QB(i), Ep8_errorFun_QB(i), Ep2_errorFun_Q(i), Ep4_errorFun_Q(i), Ep6_errorFun_Q(i), Ep8_errorFun_Q(i));
    end
    fclose(fileID);
    type(fileName);
end

end

function [Ex] = EgammaMinusNumerical(xd, f, g, h, eta)

syms x;

a = -1/2 * (g{1} + kronPolyEval(g(2:end), x)) ^ 2;
b = -kronPolyEval(f, x);
c = eta / 2 * kronPolyEval(h, x) ^ 2;

dEx2 = (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a);
dEx1 = (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a);

i = 1; Ex = zeros(1, length(xd));
for xi = xd
    if xi < 0
        Ex(i) = vpaintegral(dEx1, x, 0, xi);
    else
        Ex(i) = vpaintegral(dEx2, x, 0, xi);
    end
    i = i + 1;
end

end
