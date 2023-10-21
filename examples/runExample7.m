function [v, w] = runExample7()
%runExample7 Runs 3D aircraft stall model
%   Usage:  [v, w] = runExample7()
%
%   runExample7() runs the aircraft stall stabilization model from
%   Garrard 1977 [1].
%
%   Outputs:
%       v,w - coefficients of the past and future energy function
%             approximations, respectively.
%
%   Reference: [1] W. L. Garrard and J. M. Jordan, “Design of nonlinear
%               automatic flight control systems,” Automatica, vol. 13,
%               no. 5, pp. 497–505, Sep. 1977,
%               doi: 10.1016/0005-1098(77)90070-x
%
%   Part of the NLbalancing repository.
%%

[f, g, h] = getSystem7();
% g = g(1);
fprintf('Running Example 7\n')

options = odeset('Events', @myEvent);

%% Define original dynamics and controllers
F = @(x) kronPolyEval(f, x);

% This is the Almubarak reduced dynamics
% g = g{1};
% G = @(x) (g);
% U2 = @(x) [0;0;0]; U3 = [0;0;0];

% This is the reduced dynamics retaining nonlinear g(x)
G = @(x) (g{1} + kronPolyEval(g(2:end), x));
U2 = @(x) [0; 0; 0]; U3 = [0; 0; 0];

% This is technically the full non-affine model that Garrard presented
% G = @(x) (g{1} + kronPolyEval(g(2:end), x));
% U2 = @(x) [.47*x(1);0;46]; U3 = [.63;0;61.4];

%% Plot x1 and control u
Q = {0, 0.25, 0, 0}; R = 1; % values should be between -\infty and 1.
fprintf('Simulating for eta=%g (gamma=%g)\n', R, 1 / sqrt(1 - R))

%  Compute the polynomial approximations to the future energy function
degree = 8;
[w] = pqr(f, g, Q, R, degree);

tspan = [0, 12];

for alpha0 = [35, 30, 27, 25]
    X0 = [pi / 180 * alpha0; 0; 0];
    figure('Position', [415 47 793.3333 347.3333]);
    subplot(1, 2, 1); hold on;
    subplot(1, 2, 2); hold on;

    legendEntries = {};
    Ts = {}; X1s = {}; X2s = {}; X3s = {}; Us = {};
    for d = 2:2:degree
        u = @(x) (- R * G(x).' * kronPolyDerivEval(w(1:d), x).' / 2);
        [t, X] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0, options);

        subplot(1, 2, 1)
        plot(t, X(:, 1) / (pi / 180))
        legendEntries{end + 1} = sprintf('order-%i controller ', d - 1);

        subplot(1, 2, 2)
        us = []; Js = [];
        for i = 1:length(X)
            us(end + 1) = u(X(i, :).');
            Js(end + 1) = 0.25 * X(i, :) * X(i, :).' + u(X(i, :).') .^ 2;
        end
        valueFun_true = trapz(t, Js) / 2;
        fprintf('%i  &  %f   \n', d, valueFun_true)
        plot(t, us / (pi / 180))

        subplot(1, 2, 1)
        legend(legendEntries)
        ylim([0 alpha0 + 5])
        xlim(tspan)

        subplot(1, 2, 2)
        legend(legendEntries)
        ylim([-20 40])
        xlim(tspan)

        Ts{end + 1} = t; X1s{end + 1} = X(:, 1) / (pi / 180); Us{end + 1} = us / (pi / 180);
        X2s{end + 1} = X(:, 2) / (pi / 180); X3s{end + 1} = X(:, 3) / (pi / 180);

    end

    if exportData
        fprintf('Writing data to plots/example7_alpha%i_x1.dat \n', alpha0)
        fileID = fopen(sprintf('plots/example7_alpha%i_x1.dat', alpha0), 'w');
        fprintf(fileID, '# Figure X-a Data\n');
        fprintf(fileID, '# aircraft stall angle, angle of attack data\n');

        % Calculate the maximum number of points in any line
        max_points = max(cellfun(@numel, Ts));

        % Determine the number of lines (sets of points)
        num_lines = length(Ts);

        % Write the header
        fprintf(fileID, '       t01     &      x01      & ');

        % Write the rest of the header
        for i = 2:num_lines - 1
            fprintf(fileID, '      t%02d     &      x%02d      & ', i, i);
        end
        fprintf(fileID, '      t%02d     &      x%02d      \n ', i + 1, i + 1);

        % Iterate over the number of points
        for j = 1:max_points
            % Iterate over each line
            for i = 1:num_lines
                t_line = Ts{i};
                x_line = X1s{i};
                if j <= numel(t_line) && abs(x_line(j)) < 90 % If this line has j or more points
                    fprintf(fileID, '%+1.6e & %+1.6e', t_line(j), x_line(j));
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
        fclose(fileID); fprintf('Writing data to plots/example7_alpha%i_x1.dat \n', alpha0)
        fileID = fopen(sprintf('plots/example7_alpha%i_x1.dat', alpha0), 'w');
        fprintf(fileID, '# Figure X-a Data\n');
        fprintf(fileID, '# aircraft stall angle, angle of attack data\n');

        % Calculate the maximum number of points in any line
        max_points = max(cellfun(@numel, Ts));

        % Determine the number of lines (sets of points)
        num_lines = length(Ts);

        % Write the header
        fprintf(fileID, '       t01     &      x01      & ');

        % Write the rest of the header
        for i = 2:num_lines - 1
            fprintf(fileID, '      t%02d     &      x%02d      & ', i, i);
        end
        fprintf(fileID, '      t%02d     &      x%02d      \n ', i + 1, i + 1);

        % Iterate over the number of points
        for j = 1:max_points
            % Iterate over each line
            for i = 1:num_lines
                t_line = Ts{i};
                x_line = X1s{i};
                if j <= numel(t_line) && abs(x_line(j)) < 90 % If this line has j or more points
                    fprintf(fileID, '%+1.6e & %+1.6e', t_line(j), x_line(j));
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

        fprintf('Writing data to plots/example7_alpha%i_x2.dat \n', alpha0)
        fileID = fopen(sprintf('plots/example7_alpha%i_x2.dat', alpha0), 'w');
        fprintf(fileID, '# Figure X-ab Data\n');
        fprintf(fileID, '# aircraft stall angle, plane rotation data\n');

        % Calculate the maximum number of points in any line
        max_points = max(cellfun(@numel, Ts));

        % Determine the number of lines (sets of points)
        num_lines = length(Ts);

        % Write the header
        fprintf(fileID, '       t01     &      x01      & ');

        % Write the rest of the header
        for i = 2:num_lines - 1
            fprintf(fileID, '      t%02d     &      x%02d      & ', i, i);
        end
        fprintf(fileID, '      t%02d     &      x%02d      \n ', i + 1, i + 1);

        % Iterate over the number of points
        for j = 1:max_points
            % Iterate over each line
            for i = 1:num_lines
                t_line = Ts{i};
                x_line = X2s{i};
                if j <= numel(t_line) && abs(x_line(j)) < 90 % If this line has j or more points
                    fprintf(fileID, '%+1.6e & %+1.6e', t_line(j), x_line(j));
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

        fprintf('Writing data to plots/example7_alpha%i_x3.dat \n', alpha0)
        fileID = fopen(sprintf('plots/example7_alpha%i_x3.dat', alpha0), 'w');
        fprintf(fileID, '# Figure X-ac Data\n');
        fprintf(fileID, '# aircraft stall angle, plane rotation rate\n');

        % Calculate the maximum number of points in any line
        max_points = max(cellfun(@numel, Ts));

        % Determine the number of lines (sets of points)
        num_lines = length(Ts);

        % Write the header
        fprintf(fileID, '       t01     &      x01      & ');

        % Write the rest of the header
        for i = 2:num_lines - 1
            fprintf(fileID, '      t%02d     &      x%02d      & ', i, i);
        end
        fprintf(fileID, '      t%02d     &      x%02d      \n ', i + 1, i + 1);

        % Iterate over the number of points
        for j = 1:max_points
            % Iterate over each line
            for i = 1:num_lines
                t_line = Ts{i};
                x_line = X3s{i};
                if j <= numel(t_line) && abs(x_line(j)) < 90 % If this line has j or more points
                    fprintf(fileID, '%+1.6e & %+1.6e', t_line(j), x_line(j));
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

        fprintf('Writing data to plots/example7_alpha%i_u.dat \n', alpha0)
        fileID = fopen(sprintf('plots/example7_alpha%i_u.dat', alpha0), 'w');
        fprintf(fileID, '# Figure X-b Data\n');
        fprintf(fileID, '# aircraft stall angle, control input data\n');

        % Calculate the maximum number of points in any line
        max_points = max(cellfun(@numel, Ts));

        % Determine the number of lines (sets of points)
        num_lines = length(Ts);

        % Write the header
        fprintf(fileID, '       t01     &      x01      & ');

        % Write the rest of the header
        for i = 2:num_lines - 1
            fprintf(fileID, '      t%02d     &      x%02d      & ', i, i);
        end
        fprintf(fileID, '      t%02d     &      x%02d      \n ', i + 1, i + 1);

        % Iterate over the number of points
        for j = 1:max_points
            % Iterate over each line
            for i = 1:num_lines
                t_line = Ts{i};
                x_line = Us{i};
                if j <= numel(t_line) % If this line has j or more points
                    fprintf(fileID, '%+1.6e & %+1.6e', t_line(j), x_line(j));
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

        fprintf('Writing data to plots/example7_alpha%i_du.dat \n', alpha0)
        fileID = fopen(sprintf('plots/example7_alpha%i_du.dat', alpha0), 'w');
        fprintf(fileID, '# Figure X-b Data\n');
        fprintf(fileID, '# aircraft stall angle, control input derivative data\n');

        % Calculate the maximum number of points in any line
        max_points = max(cellfun(@numel, Ts)) - 1;

        % Determine the number of lines (sets of points)
        num_lines = length(Ts);

        % Write the header
        fprintf(fileID, '       t01     &      x01      & ');

        % Write the rest of the header
        for i = 2:num_lines - 1
            fprintf(fileID, '      t%02d     &      x%02d      & ', i, i);
        end
        fprintf(fileID, '      t%02d     &      x%02d      \n ', i + 1, i + 1);

        % Iterate over the number of points
        for j = 1:max_points
            % Iterate over each line
            for i = 1:num_lines
                t_line = Ts{i};
                x_line = diff(Us{i}) ./ diff(Ts{i}.');
                if j <= numel(t_line) - 1 % If this line has j or more points
                    fprintf(fileID, '%+1.6e & %+1.6e', t_line(j), x_line(j));
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

end

end

function [value, isterminal, direction] = myEvent(T, Y)
value = (abs(Y(1)) > 90 * pi / 180);
isterminal = 1; % Stop the integration
direction = 0;
end
