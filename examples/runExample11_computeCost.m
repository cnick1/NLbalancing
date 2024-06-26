function [v, w] = runExample11_computeCost(exportPlotData, nFterms, degree, varargin)
%runExample11_computeCost Runs the 2D example to plot energy functions as contour plots
%
%   Usage:  [v,w] = runExample11_computeCost(degree,plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel)
%
%   runExample11_computeCost() runs the default case of a quadratic model from [1] which
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
if nargin < 3
    if nargin < 2
        if nargin < 1
            exportPlotData = false;
        end
        nFterms = 5;
    end
    degree = 8;
end

Q = 0;
R = 1; Rinv = inv(R);

if nFterms == 1
    nFterms = 2; % Note F2 is zero, this is just to be able to compute a controller and ignore the error if F2 doesn't exist
end

%% Get model and compute energy functions
scale = .1767; scaling = 1 / sqrt(scale); % For plot and initial condition scaling, hardcoded

m = 1; L = 10; %56.5962*scale;
gravity = 9.81;
[f, g, h] = getSystem11(nFterms, m, L);
fprintf('Running Example 11\n')

%% Closed-loop phase portraits
% Define the range of initial conditions

y0 = [sqrt(6 * gravity / L)]; x0 = [-pi];
% y0 = [sqrt(3*gravity/L)]; x0 = [-pi/2];
% y0 = [sqrt(3*gravity/L*(1-sqrt(2)/2))]; x0 = [-pi/4];

tspan = [0 10];

%  Compute the polynomial approximations to the past future energy function
% [w] = ppr(f, g, h2q(h), eta, degree);
options.verbose = true;
[w] = ppr(f, g, Q, R, degree, options);

syms x1 x2
vpa(-Rinv * g{1}.' * kronPolyDerivEval(w, [x1; x2]).') / 2

options = odeset('Events', @myEvent);

xsNL = {}; ysNL = {};
[t, y] = ode45(@(t, y) [y(2); 3 * gravity / (2 * L) * sin(y(1))] - Rinv * g{1} * g{1}.' * (0.5 * kronPolyDerivEval(w, y).'), tspan, [x0; y0], options);
u = zeros(length(y), 1);
for i = 1:length(y)
    u(i) =- Rinv * g{1}.' * (0.5 * kronPolyDerivEval(w, y(i, :).').');
end

% Create a figure and set up subplots
figure
subplot(1, 2, 1); hold on;
plot(y(:, 1), y(:, 2), 'r'); xsNL{end + 1} = y(:, 1); ysNL{end + 1} = y(:, 2);
subplot(1, 2, 2); hold on;
plot(t, u, 'b');

l2_norm_valueFun = (0.5 * kronPolyEval(w, [x0; y0]));
l2_norm_true = trapz(t, u .^ 2) / 2;

fprintf('%i  &  %f  &  %f   \n', degree, l2_norm_valueFun, l2_norm_true)

% compute L2 norm of the control

% Set up the plot
xlim([-pi pi]); ylim([-2 * scaling 2 * scaling]); xlabel('x'); title('closed loop pendulum');

if false %exportPlotData
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
