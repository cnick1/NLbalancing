function [v, w] = runExample1(numGTermsModel, numGTermsApprox, exportPlotData, varargin)
%runExample1 Runs 1D ODE example to compare computed and analytical energy functions
%   Usage:  [v, w] = runExample1(numGTermsModel, numGTermsApprox, exportPlotData)
%
%   runExample1() runs the default case of a quadratic model from [1].
%
%   Inputs:
%       numGTermsModel     - Number of terms in the full order model
%       numGTermsApprox    - Number of terms assumed when computing energy functions
%       exportPlotData     - Boolean variable to determine if plots/data are exported
%
%   Outputs:
%       v,w             are coefficients of the past and future energy
%                       function approximations, respectively.
%
%   Reference: [1] Nonlinear Balanced Truncation Model Reduction:
%        Part 1-Computing Energy Functions, by Kramer, Gugercin, and Borggaard.
%        arXiv:2209.07645.
%
%   Part of the NLbalancing repository.
%%

if nargin < 3
    exportPlotData = false;
    if nargin < 2
        if nargin < 1
            numGTermsModel = 1;
        end
        numGTermsApprox = numGTermsModel;
    end
end

[A, B, C, N, f, g, h] = getSystem1();
g(numGTermsModel + 1:end) = deal({0}); % Adjust FOM to be Quadratic, QB, etc.
fprintf('Running Example 1\n')

eta = 0.5; % values should be between -\infty and 1.
% eta=0.5 corresponds to gamma= sqrt(2)
% since eta = 1 - 1/gamma^2;

fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

%  Compute the polynomial approximations to the future energy function
d = 8;
[w, RES] = approxFutureEnergy(f, N, g(1:numGTermsApprox), C, eta, d, true);

w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};

x = linspace(-6, 6, 250);

Ef2 = 0.5 * w2 * x .^ 2;
Ef3 = Ef2 + 0.5 * w3 * x .^ 3;
Ef4 = Ef3 + 0.5 * w4 * x .^ 4;
Ef5 = Ef4 + 0.5 * w5 * x .^ 5;
Ef6 = Ef5 + 0.5 * w6 * x .^ 6;
Ef7 = Ef6 + 0.5 * w7 * x .^ 7;
Ef8 = Ef7 + 0.5 * w8 * x .^ 8;

%   Efd = Ef2;
%   for idx = 3:length(w)
%     Efd = Efd + 0.5*w{idx}*x.^idx;
%   end

%  Compute the analytical solution for comparison

EPlusAnalytic = EgammaPlusNumerical(x,f,g,h, eta);

figure
plot(x(1:10:end), EPlusAnalytic(1:10:end), '+', ...
    x, Ef2, ...
    x, Ef4, ...
    x, Ef6, ...
    x, Ef8, ...
    'LineWidth', 2)
legend('analytic', ...
    'degree 2', ...
    'degree 4', ...
    'degree 6', ...
    'degree 8', ...
    'Location', 'northwest')
xlabel('$x$', ...
    'interpreter', 'latex', ...
    'FontSize', 20, ...
    'fontweight', 'bold')
ylabel('$\mathcal{E}_\gamma^+$', ...
    'interpreter', 'latex', ...
    'FontSize', 20, ...
    'fontweight', 'bold')

% xlim([-3, 3])
% ylim([0, 3])

if exportPlotData
    %  Save data to generate tikz plots for the paper
    fid = fopen('plots/ex1_future_a.txt', 'w');
    fprintf(fid, '%g %g\n', [x; EPlusAnalytic]);
    fclose(fid);

    fid = fopen('plots/ex1_future_2.txt', 'w');
    fprintf(fid, '%g %g\n', [x; Ef2]);
    fclose(fid);

    fid = fopen('plots/ex1_future_4.txt', 'w');
    fprintf(fid, '%g %g\n', [x; Ef4]);
    fclose(fid);

    fid = fopen('plots/ex1_future_6.txt', 'w');
    fprintf(fid, '%g %g\n', [x; Ef6]);
    fclose(fid);

    fid = fopen('plots/ex1_future_8.txt', 'w');
    fprintf(fid, '%g %g\n', [x; Ef8]);
    fclose(fid);
end
%   %  Compute the polynomial approximations to the past energy function
[v] = approxPastEnergy(f, N, g(1:numGTermsApprox), C, eta, d);
v2 = v{2}; v3 = v{3}; v4 = v{4}; v5 = v{5}; v6 = v{6}; v7 = v{7}; v8 = v{8};

% x = linspace(-6,6,121);

Ep2 = 0.5 * v2 * x .^ 2;
Ep3 = Ep2 + 0.5 * v3 * x .^ 3;
Ep4 = Ep3 + 0.5 * v4 * x .^ 4;
Ep5 = Ep4 + 0.5 * v5 * x .^ 5;
Ep6 = Ep5 + 0.5 * v6 * x .^ 6;
Ep7 = Ep6 + 0.5 * v7 * x .^ 7;
Ep8 = Ep7 + 0.5 * v8 * x .^ 8;

%  Compute the analytical solution for comparison
EMinusAnalytic = EgammaMinusNumerical(x, f, g, h, eta);

figure
plot(x(1:10:end), EMinusAnalytic(1:10:end), '+', ...
    x, Ep2, ...
    x, Ep4, ...
    x, Ep6, ...
    x, Ep8, ...
    'LineWidth', 2)
legend('analytic', ...
    'degree 2', ...
    'degree 4', ...
    'degree 6', ...
    'degree 8', ...
    'Location', 'northeast')
xlabel('$x$', ...
    'interpreter', 'latex', ...
    'FontSize', 20, ...
    'fontweight', 'bold')
ylabel('$\mathcal{E}_\gamma^-$', ...
    'interpreter', 'latex', ...
    'FontSize', 20, ...
    'fontweight', 'bold')

if exportPlotData
    %   fid = fopen('plots/ex1_errorpe_2.txt','w');
    %   fprintf(fid,'%g %g\n',[x;abs(Ep2-EMinusAnalytic)]);
    %   fclose(fid);
    %
    %   fid = fopen('plots/ex1_errorpe_4.txt','w');
    %   fprintf(fid,'%g %g\n',[x;abs(Ep4-EMinusAnalytic)]);
    %   fclose(fid);
    %
    %   fid = fopen('plots/ex1_errorpe_6.txt','w');
    %   fprintf(fid,'%g %g\n',[x;abs(Ep6-EMinusAnalytic)]);
    %   fclose(fid);
    %
    %   fid = fopen('plots/ex1_errorpe_8.txt','w');
    %   fprintf(fid,'%g %g\n',[x;abs(Ep8-EMinusAnalytic)]);
    %   fclose(fid);
    %
    %   %
    fid = fopen('plots/ex1_past_a.txt', 'w');
    fprintf(fid, '%g %g\n', [x; EMinusAnalytic]);
    fclose(fid);

    fid = fopen('plots/ex1_past_2.txt', 'w');
    fprintf(fid, '%g %g\n', [x; Ep2]);
    fclose(fid);

    fid = fopen('plots/ex1_past_4.txt', 'w');
    fprintf(fid, '%g %g\n', [x; Ep4]);
    fclose(fid);

    fid = fopen('plots/ex1_past_6.txt', 'w');
    fprintf(fid, '%g %g\n', [x; Ep6]);
    fclose(fid);

    fid = fopen('plots/ex1_past_8.txt', 'w');
    fprintf(fid, '%g %g\n', [x; Ep8]);
    fclose(fid);
    %
    %
    %   save('Ex1_RawData.mat','x',...
    %        'EPlusAnalytic', 'Ef2','Ef3','Ef4','Ef5','Ef6','Ef7','Ef8',...
    %        'EMinusAnalytic','Ep2','Ep3','Ep4','Ep5','Ep6','Ep7','Ep8')
end
end

function [Ex] = EgammaPlusNumerical(xd, f, g, h, eta)

syms x;

a = -eta/2 * (g{1} + kronPolyEval(g(2:end),x)) ^ 2;
b = kronPolyEval(f,x);
c = 1/2 * kronPolyEval(h,x) ^ 2;

dEx2 = (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a);
dEx1 = (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a);

% Ex1 = int(dEx1, x);
% Ex2 = int(dEx2, x);

% y = piecewise(x < 0, Ex1, x >= 0, Ex2);

i = 1; Ex = zeros(1, length(xd));
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

function [Ex] = EgammaMinusNumerical(xd, f, g, h, eta)

syms x;

a = -1/2 * (g{1} + kronPolyEval(g(2:end),x)) ^ 2;
b = -kronPolyEval(f,x);
c = eta/2 * kronPolyEval(h,x) ^ 2;

dEx2 = (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a);
dEx1 = (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a);

% Ex1 = int(dEx1, x);
% Ex2 = int(dEx2, x);

% y = piecewise(x < 0, Ex1, x >= 0, Ex2);

i = 1; Ex = zeros(1, length(xd));
for xi = xd
    if xi < 0
        Ex(i) = vpaintegral(dEx1, x, 0, xi);
    else
        Ex(i) = vpaintegral(dEx2, x, 0, xi);
    end
    i = i + 1;
end

% hold on
% plot(xd(1:10:end),Ex(1:10:end),'+', 'LineWidth', 1)

end


