function [v, w] = runExample11(exportPlotData, nFterms, degree, eta, varargin)
%runExample11 Runs the 2D example to plot energy functions as contour plots
%
%   Usage:  [v,w] = runExample11(degree,plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel)
%
%   runExample11() runs the default case of a quadratic model from [1] which
%                 is based on a model from [2].
%
%   Inputs:
%       degree          is the degree of energy function approximations
%       exportPlotData   Boolean variable to determine if plots/data are exported
%
%
%   Outputs:
%       v,w              are coefficients of the past and future energy
%                        function approximations, respectively.
%
%   The value of eta is set below.
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%% Process inputs
if nargin < 4
    if nargin < 3
        if nargin < 2
            if nargin < 1
                exportPlotData = false;
            end
            nFterms = 5;
        end
        degree = 8;
    end
    % Compute energy functions
    eta = 1; % values should be between -\infty and 1.
    % eta=1 is HJB/closed-loop balancing, 0 is open loop.
end

if nFterms == 1
    nFterms = 2; % Note F2 is zero, this is just to be able to compute a controller and ignore the error if F2 doesn't exist
end

%% Get model and compute energy functions
scale = .1767; scaling = 1 / sqrt(scale); % For plot and initial condition scaling, hardcoded

m = 1; L = 10; %56.5962*scale;
gravity = 9.81;
[f, g, h] = getSystem11(nFterms, m, L);
fprintf('Running Example 11\n')

%% Open loop phase portraits
% Define the range of initial conditions
y0 = [-2, -1.7, -1.5, -1.3, -1.1, - .9, - .8, - .6, - .4, - .2, - .05, - .01, 0., 0.01, 0.05, 0.2, 0.4, 0.6, 0.8, .95, 1.1, 1.3, 1.5, 1.7, 2] * scaling;
x0 = [zeros(1, 11) + pi, zeros(1, 3), zeros(1, 11) - pi];

y0 = [y0, sqrt(6 * gravity / L)];
x0 = [x0, -pi];

tspan = [0 7];

% Create a figure and set up subplots
figure('Position', [600 60 860 240]); subplot(1, 3, 1); hold on;

xsNL = {}; ysNL = {};
% Loop over the initial conditions and solve the ODE
for i = 1:length(y0) - 1
    [t, y] = ode45(@(t, y) [y(2); 3 * gravity / (2 * L) * sin(y(1))], tspan, [x0(i); y0(i)]);
    plot(y(:, 1), y(:, 2), 'r'); xsNL{end + 1} = y(:, 1); ysNL{end + 1} = y(:, 2);
    [t, y] = ode45(@(t, y) [y(2); 3 * gravity / (2 * L) * sin(y(1))], -tspan, [x0(i); y0(i)]);
    plot(y(:, 1), y(:, 2), 'r'); xsNL{end + 1} = y(:, 1); ysNL{end + 1} = y(:, 2);
end
i = length(y0);
[t, y] = ode23(@(t, y) [y(2); 3 * gravity / (2 * L) * sin(y(1))], tspan, [x0(i); y0(i)]);
plot(y(:, 1), y(:, 2), 'b'); xsNL{end + 1} = y(:, 1); ysNL{end + 1} = y(:, 2);
[t, y] = ode23(@(t, y) [y(2); 3 * gravity / (2 * L) * sin(y(1))], -tspan, [x0(i); y0(i)]);
plot(y(:, 1), y(:, 2), 'b'); xsNL{end + 1} = y(:, 1); ysNL{end + 1} = y(:, 2);

% Set up the plot
xlim([-pi pi]); ylim([-2 * scaling 2 * scaling]); xlabel('x'); title('pendulum');

% Set up the second subplot
subplot(1, 3, 2); hold on;

xsPoly = {}; ysPoly = {};
% Loop over the initial conditions and solve the ODE
for i = 1:length(y0) - 1
    [t, y] = ode45(@(t, x) kronPolyEval(f, x), tspan, [x0(i); y0(i)]);
    plot(y(:, 1), y(:, 2), 'r'); xsPoly{end + 1} = y(:, 1); ysPoly{end + 1} = y(:, 2);
    [t, y] = ode45(@(t, x) kronPolyEval(f, x), -tspan, [x0(i); y0(i)]);
    plot(y(:, 1), y(:, 2), 'r'); xsPoly{end + 1} = y(:, 1); ysPoly{end + 1} = y(:, 2);
end

i = length(y0);
[t, y] = ode45(@(t, x) kronPolyEval(f, x), tspan, [x0(i); y0(i)]);
plot(y(:, 1), y(:, 2), 'b'); xsPoly{end + 1} = y(:, 1); ysPoly{end + 1} = y(:, 2);
[t, y] = ode45(@(t, x) kronPolyEval(f, x), -tspan, [x0(i); y0(i)]);
plot(y(:, 1), y(:, 2), 'b'); xsPoly{end + 1} = y(:, 1); ysPoly{end + 1} = y(:, 2);

% Set up the plot
xlim([-pi pi]); ylim([-2 * scaling 2 * scaling]); xlabel('x'); title('polynomial approx.');

if exportPlotData
    xs = xsNL; ys = ysNL;
    fprintf('Writing data to plots/example11_openLoopPhasePortraits_nonlinear.dat \n')
    fileID = fopen('plots/example11_openLoopPhasePortraits_nonlinear.dat', 'w');
    fprintf(fileID, '# Figure X-a Data\n');
    fprintf(fileID, '# pendulum open loop phase portrait trajectory data\n');

    % Calculate the maximum number of points in any line
    max_points = max(cellfun(@numel, xs));

    % Determine the number of lines (sets of points)
    num_lines = length(xs);

    % Write the header
    fprintf(fileID, '       x01     &      y01      & ');

    % Write the rest of the header
    for i = 2:num_lines - 1
        fprintf(fileID, '      x%02d     &      y%02d      & ', i, i);
    end
    fprintf(fileID, '      x%02d     &      y%02d      \n ', i + 1, i + 1);

    % Iterate over the number of points
    for j = 1:max_points
        % Iterate over each line
        for i = 1:num_lines
            x_line = xs{i};
            y_line = ys{i};
            if j <= numel(x_line) && abs(x_line(j)) < 4.5 % If this line has j or more points
                fprintf(fileID, '%+1.6e & %+1.6e', x_line(j), y_line(j));
                % Add '&' delimiter unless it's the last line
                if i < num_lines
                    fprintf(fileID, ' & ');
                else
                    fprintf(fileID, ' \n ');
                end
            else % this line doesn't have a jth point
                fprintf(fileID, '              &              ');
                % Add '&' delimiter unless it's the last line
                if i < num_lines
                    fprintf(fileID, ' & ');
                else
                    fprintf(fileID, ' \n ');
                end
            end
        end
    end
    % Close the data file
    fclose(fileID);

    xs = xsPoly; ys = ysPoly;
    fprintf('Writing data to plots/example11_openLoopPhasePortraits_polynomial%i.dat \n', nFterms)
    fileID = fopen(sprintf('plots/example11_openLoopPhasePortraits_polynomial%i.dat', nFterms), 'w');
    fprintf(fileID, '# Figure X-b Data\n');
    fprintf(fileID, '# pendulum open loop phase portrait trajectory data, degree %i approximation\n', nFterms);

    % Calculate the maximum number of points in any line
    max_points = max(cellfun(@numel, xs));

    % Determine the number of lines (sets of points)
    num_lines = length(xs);

    % Write the header
    fprintf(fileID, '       x01     &      y01      & ');

    % Write the rest of the header
    for i = 2:num_lines - 1
        fprintf(fileID, '      x%02d     &      y%02d      & ', i, i);
    end
    fprintf(fileID, '      x%02d     &      y%02d      \n ', i + 1, i + 1);

    % Iterate over the number of points
    for j = 1:max_points
        % Iterate over each line
        for i = 1:num_lines
            x_line = xs{i};
            y_line = ys{i};
            if j <= numel(x_line) && abs(x_line(j)) < 4.5 % If this line has j or more points
                fprintf(fileID, '%+1.6e & %+1.6e', x_line(j), y_line(j));
                % Add '&' delimiter unless it's the last line
                if i < num_lines
                    fprintf(fileID, ' & ');
                else
                    fprintf(fileID, ' \n ');
                end
            else % this line doesn't have a jth point
                fprintf(fileID, '              &              ');
                % Add '&' delimiter unless it's the last line
                if i < num_lines
                    fprintf(fileID, ' & ');
                else
                    fprintf(fileID, ' \n ');
                end
            end
        end
    end
    % Close the data file
    fclose(fileID);
end

%% Closed loop phase portraits
% Define the range of initial conditions
% x0 = [-2, -1.7, -1.5, -1.3, -1.1, - .9, - .8, - .6, - .4, - .2, - .05, - .01, 0., 0.01, 0.05, 0.2, 0.4, 0.6, 0.8, .95, 1.1, 1.3, 1.5, 1.7, 2] * scaling;
% y0 = [zeros(1, 11), zeros(1, 3), zeros(1, 11)];

fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

%  Compute the polynomial approximations to the past future energy function
% [v] = approxPastEnergy(f, N, g, h, eta, degree, true);
% [w] = pqr(f, g, h2q(h), eta, degree, true);
[w] = pqr(f, g, 0, eta, degree, true);

% Create a figure and set up subplots
subplot(1, 3, 3); hold on;

options = odeset('Events', @myEvent);

xsNL = {}; ysNL = {};
% Loop over the initial conditions and solve the ODE
for i = 1:length(y0)
    [t, y] = ode45(@(t, y) [y(2); 3 * gravity / (2 * L) * sin(y(1))] - eta * g{1} * g{1}.' * (0.5 * kronPolyDerivEval(w, y).'), tspan, [x0(i); y0(i)], options);
    plot(y(:, 1), y(:, 2), 'r'); xsNL{end + 1} = y(:, 1); ysNL{end + 1} = y(:, 2);
    [t, y] = ode45(@(t, y) [y(2); 3 * gravity / (2 * L) * sin(y(1))] - eta * g{1} * g{1}.' * (0.5 * kronPolyDerivEval(w, y).'), -tspan, [x0(i); y0(i)], options);
    plot(y(:, 1), y(:, 2), 'r'); xsNL{end + 1} = y(:, 1); ysNL{end + 1} = y(:, 2);
end

i = length(y0);
[t, y] = ode45(@(t, y) [y(2); 3 * gravity / (2 * L) * sin(y(1))] - eta * g{1} * g{1}.' * (0.5 * kronPolyDerivEval(w, y).'), tspan, [x0(i); y0(i)], options);
plot(y(:, 1), y(:, 2), 'b'); xsNL{end + 1} = y(:, 1); ysNL{end + 1} = y(:, 2);
[t, y] = ode45(@(t, y) [y(2); 3 * gravity / (2 * L) * sin(y(1))] - eta * g{1} * g{1}.' * (0.5 * kronPolyDerivEval(w, y).'), -tspan, [x0(i); y0(i)], options);
plot(y(:, 1), y(:, 2), 'b'); xsNL{end + 1} = y(:, 1); ysNL{end + 1} = y(:, 2);

% Set up the plot
xlim([-pi pi]); ylim([-2 * scaling 2 * scaling]); xlabel('x'); title('closed loop pendulum');

if exportPlotData
    xs = xsNL; ys = ysNL;
    fprintf('Writing data to plots/example11_closedLoopPhasePortraits_d%i_polynomial%i.dat \n', degree, nFterms)
    fileID = fopen(sprintf('plots/example11_closedLoopPhasePortraits_d%i_polynomial%i.dat', degree, nFterms), 'w');
    fprintf(fileID, '# Figure X-a Data\n');
    fprintf(fileID, '# pendulum closed loop phase portrait trajectory data\n');

    % Calculate the maximum number of points in any line
    max_points = max(cellfun(@numel, xs));

    % Determine the number of lines (sets of points)
    num_lines = length(xs);

    % Write the header
    fprintf(fileID, '       x01     &      y01      & ');

    % Write the rest of the header
    for i = 2:num_lines - 1
        fprintf(fileID, '      x%02d     &      y%02d      & ', i, i);
    end
    fprintf(fileID, '      x%02d     &      y%02d      \n ', i + 1, i + 1);

    % Iterate over the number of points
    for j = 1:max_points
        if exist('count', 'var') == 1 && count == 50
            break
        else
            count = 0;
        end
        % Iterate over each line
        for i = 1:num_lines
            x_line = xs{i};
            y_line = ys{i};
            if j <= numel(x_line) && abs(x_line(j)) < 4.5 % If this line has j or more points
                fprintf(fileID, '%+1.6e & %+1.6e', x_line(j), y_line(j));
                % Add '&' delimiter unless it's the last line
                if i < num_lines
                    fprintf(fileID, ' & ');
                else
                    fprintf(fileID, ' \n ');
                end
            else % this line doesn't have a jth point
                count = count + 1;
                fprintf(fileID, '              &              ');
                % Add '&' delimiter unless it's the last line
                if i < num_lines
                    fprintf(fileID, ' & ');
                else
                    fprintf(fileID, ' \n ');
                end
            end
        end
    end
    % Close the data file
    fclose(fileID);

end

end

function [value, isterminal, direction] = myEvent(T, Y)
value = max(abs(Y)) > 15; % check if any element of Y is greater than 1e6
isterminal = 1; % stop the integration if value is true
direction = 0; % direction doesn't matter
end
