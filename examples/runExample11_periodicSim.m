function [v, w] = runExample11_periodicSim(exportPlotData, nFterms, degree, eta, varargin)
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
            nFterms = 9;
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

m = 1; L = 10; gravity = 9.81;
[f, g, h] = getSystem11(nFterms, m, L);

%% Closed loop phase portraits
% Define the range of initial conditions
tspan = [0 7];

x0 = [-2];
y0 = [-4];

[w] = pqr(f, g, h, 1 / eta, degree, true);

% Create a figure and set up subplots
figure; hold on;

% Loop over the initial conditions and solve the ODE
for i = 1:length(y0)
    [t, y] = ode45(@(t, y) [y(2); 3 * gravity / (2 * L) * sin(y(1))], tspan, [x0(i); y0(i)]);
    y(:, 1) = mod(y(:, 1) + pi, 2 * pi) - pi;
    [X, Y] = parseTrajectories(y);
    plot(X, Y, 'k', 'LineWidth', 2);
    [t, y] = ode45(@(t, y) kronPolyEval(f, y), tspan, [x0(i); y0(i)]);
    plot(y(:, 1), y(:, 2), 'r');
    %     y(:,1) = mod(y(:,1) + pi, 2*pi) - pi;
    %     [X, Y] = parseTrajectories(y);
    plot(X, Y, 'r-.', 'LineWidth', 2);
    [t, y] = ode45(@(t, y) kronPolyEvalMod(f, y), tspan, [x0(i); y0(i)]); y(:, 1) = mod(y(:, 1) + pi, 2 * pi) - pi;
    [X, Y] = parseTrajectories(y);
    plot(X, Y, 'g--', 'LineWidth', 2);
end

% Set up the plot
xlim([-pi pi]); ylim([-2 * scaling 2 * scaling]); xlabel('x'); title('closed loop pendulum');

end

function FofY = kronPolyEvalMod(f, y)
y(1) = mod(y(1) + pi, 2 * pi) - pi;
FofY = kronPolyEval(f, y);
end

function [X, Y] = parseTrajectories(y)

% Find the indices where the trajectory crosses the boundary
boundary_crossings = find(abs(diff(y(:, 1))) > pi);

max_length = length(y(:, 1));
X = NaN(max_length, 1); Y = X;

% Initialize plot variables
start_idx = 1;

% Loop through each segment
for i = 1:length(boundary_crossings) + 1
    if i <= length(boundary_crossings)
        end_idx = boundary_crossings(i);
    else
        end_idx = size(y, 1);
    end

    % Plot the current segment
    %     plot(y(start_idx:end_idx, 1), y(start_idx:end_idx, 2), 'b-'); % Adjust line style as needed

    % Pad the shorter vector with NaNs to match the maximum length
    X = [X, [y(start_idx:end_idx, 1); NaN(max_length - length(y(start_idx:end_idx, 1)), 1)]];
    Y = [Y, [y(start_idx:end_idx, 2); NaN(max_length - length(y(start_idx:end_idx, 2)), 1)]];

    % Update the starting index for the next segment
    start_idx = end_idx + 1;
end

end
