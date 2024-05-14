function runExample6(exportData, x0)
%runExample6 Runs the finite element beam example to demonstrate
%            convergence and scalability.
%
%   Usage:  runExample6(exportData, x0)
%
%   Inputs:
%       exportData - Boolean variable to determine if
%                         plots/data are exported
%
%   Outputs:
%       v,w        - Coefficients of the past and future energy
%                         function approximations, respectively
%
%   The value of eta is set below.
%
%   Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” 2023.
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 6\n')

if nargin < 2
    if nargin < 1
        exportData = false;
    end
    x0 = 0.01;
end

eta = 0.1; 

%% Table II: convergence & scalability wrt n (d=3)
%  Since the initial times are so short, we average nTest times

degree = 3;

% print results to command window as we run, then write to file
fprintf('\n# Table II Data\n, # finite element beam model, convergence and scalability results; d=%d \nnumElements &   n  &     n^%d       &  CPU-sec  &  E_%d^+(x_0)     \n', degree, degree, degree);

nTest = 10;
nd = []; times = []; energies = [];
numEls = [1, 2, 4, 8, 16];
for numEl = numEls
    fprintf(' %5d      &%4d  & ', numEl, 6 * numEl);

    % Compute energy functions and CPU time
    [f, g, h] = getSystem6(numEl, 2);
    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;
    
    % Evaluate energy function at x0 corresponding to nodes having linear
    % displacement but no rotation or initial velocity
    numNodes = numEl + 1;
    initialCondition = x0 / (numNodes - 1) * ...
        [[(0:numNodes - 1);
        (0:numNodes - 1);
        0 * (0:numNodes - 1)].';
        zeros(numNodes, 3)].'; % Full initial condition
    initialCondition(:, 1) = []; initialCondition(:, 1 + numNodes) = []; % Remove the first node DOFs
    initialCondition = initialCondition(:);
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);

    % Print results to command window
    fprintf(' %10.4e   & %8.2e  & %12.6e    \n', length(w{degree}), tt, wzInit);

    % log results to write to file later
    nd = [nd, length(w{degree})]; times = [times, tt]; energies = [energies, wzInit];
end

% For run-time, only run the higher cases if exporting data, and only
% run once since the run time is longer so error is less sensitive
if exportData
    nTest = 1;
    for numEl = [32, 64]%, 128] %128 runs out of ram in kroneckerLeft.m
        fprintf(' %5d      &%4d  & ', numEl, 6 * numEl);

        % Compute energy functions and CPU time
        [f, g, h] = getSystem6(numEl, 2);
        tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;

        % Evaluate energy function at x0 corresponding to nodes having linear
        % displacement but no rotation or initial velocity
        numNodes = numEl + 1;
        initialCondition = x0 / (numNodes - 1) * ...
            [[(0:numNodes - 1);
            (0:numNodes - 1);
            0 * (0:numNodes - 1)].';
            zeros(numNodes, 3)].'; % Full initial condition
        initialCondition(:, 1) = []; initialCondition(:, 1 + numNodes) = []; % Remove the first node DOFs
        initialCondition = initialCondition(:);

        wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);

        % Print results to command window
        fprintf(' %10.4e   & %8.2e  & %12.6e    \n', length(w{degree}), tt, wzInit);

        % log results to write to file later
        numEls = [numEls, numEl]; nd = [nd, length(w{degree})]; times = [times, tt]; energies = [energies, wzInit];
    end
end
if exportData
    logy = log10(times); % take the natural log of y data
    logx = log10(6 * numEls); % take the natural log of x data
    X = [ones(length(logx), 1) logx']; % create matrix of x data with a column of ones
    beta = X \ logy'; % solve for beta coefficients using linear least squares
    a = beta(1); % calculate exponent d
    d = beta(2); % calculate exponent d
    fprintf('The exponent fit gives n^%f compared with n^%d. \nWriting data to plots/example6_convergenceData_d3.dat \n', d, degree)
    fileID = fopen('plots/example6_convergenceData_d3.dat', 'w');
    fprintf(fileID, '# Table I Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)     &  exponentCoeff  &  exponentFit \n', degree, degree, degree);
    %print the header
    for i = 1:length(numEls)
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e   &  %2.2f  &  %2.2f \n', numEls(i), 6 * numEls(i), nd(i), times(i), energies(i), a, d);
    end
    fclose(fileID);
end

%% Table II: convergence & scalability wrt n (d=4)
degree = 4;

%print the header
fprintf('\n# Table II Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &   n  &     n^%d       &  CPU-sec  &  E_%d^+(x_0)     \n', degree, degree, degree);

% compute and print the results
nTest = 3;
nd = []; times = []; energies = [];
numEls = [1, 2, 4];
for numEl = numEls
    fprintf(' %5d      &%4d  & ', numEl, 6 * numEl);

    % Compute energy functions and CPU time
    [f, g, h] = getSystem6(numEl, 2);
    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;
    
    % Evaluate energy function at x0 corresponding to nodes having linear
    % displacement but no rotation or initial velocity
    numNodes = numEl + 1;
    initialCondition = x0 / (numNodes - 1) * ...
        [[(0:numNodes - 1);
        (0:numNodes - 1);
        0 * (0:numNodes - 1)].';
        zeros(numNodes, 3)].'; % Full initial condition
    initialCondition(:, 1) = []; initialCondition(:, 1 + numNodes) = []; % Remove the first node DOFs
    initialCondition = initialCondition(:);
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);

    % Print results to command window
    fprintf(' %10.4e   & %8.2e  & %12.6e    \n', length(w{degree}), tt, wzInit);

    % log results to write to file later
    nd = [nd, length(w{degree})]; times = [times, tt]; energies = [energies, wzInit];
end

% For run-time, only run the higher cases if exporting data, and only
% run once since the run time is longer so error is less sensitive
if exportData
    nTest = 1;
    for numEl = [8, 16]%, 32, 64]
        fprintf(' %5d      &%4d  & ', numEl, 6 * numEl);

        % Compute energy functions and CPU time
        [f, g, h] = getSystem6(numEl, 2);
        tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;

        % Evaluate energy function at x0 corresponding to nodes having linear
        % displacement but no rotation or initial velocity
        numNodes = numEl + 1;
        initialCondition = x0 / (numNodes - 1) * ...
            [[(0:numNodes - 1);
            (0:numNodes - 1);
            0 * (0:numNodes - 1)].';
            zeros(numNodes, 3)].'; % Full initial condition
        initialCondition(:, 1) = []; initialCondition(:, 1 + numNodes) = []; % Remove the first node DOFs
        initialCondition = initialCondition(:);

        wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);

        % Print results to command window
        fprintf(' %10.4e   & %8.2e  & %12.6e    \n', length(w{degree}), tt, wzInit);

        % log results to write to file later
        numEls = [numEls, numEl]; nd = [nd, length(w{degree})]; times = [times, tt]; energies = [energies, wzInit];
    end
end
if exportData
    logy = log10(times); % take the natural log of y data
    logx = log10(6 * numEls); % take the natural log of x data
    X = [ones(length(logx), 1) logx']; % create matrix of x data with a column of ones
    beta = X \ logy'; % solve for beta coefficients using linear least squares
    a = beta(1); % calculate exponent d
    d = beta(2); % calculate exponent d
    fprintf('The exponent fit gives n^%f compared with n^%d. \nWriting data to plots/example6_convergenceData_d4.dat \n', d, degree)
    fileID = fopen('plots/example6_convergenceData_d4.dat', 'w');
    fprintf(fileID, '# Table II Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)     &  exponentCoeff  &  exponentFit \n', degree, degree, degree);
    for i = 1:length(numEls)
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e   &  %2.2f  &  %2.2f \n', numEls(i), 6 * numEls(i), nd(i), times(i), energies(i), a, d);
    end
    fclose(fileID);
end

%% Table III: convergence & scalability wrt d (n=18)

numEl = 3;

fprintf('\n# Table III Data\n# finite element beam model, convergence and scalability results \n# numEls = %d   -->   n = %d \n   d   &  CPU-sec-2   & E_d^+(x_0)      \n', numEl, 6 * numEl);

% compute and print the results
nTest = 3;

[f, g, h] = getSystem6(numEl);

% Evaluate energy function at x0 corresponding to nodes having linear
% displacement but no rotation or initial velocity
numNodes = numEl + 1;
initialCondition = x0 / (numNodes - 1) * ...
    [[(0:numNodes - 1);
    (0:numNodes - 1);
    0 * (0:numNodes - 1)].';
    zeros(numNodes, 3)].'; % Full initial condition
initialCondition(:, 1) = []; initialCondition(:, 1 + numNodes) = []; % Remove the first node DOFs
initialCondition = initialCondition(:);

pastTimes = []; futureTimes = []; pastEnergies = []; futureEnergies = [];
degrees = 2:4;
for degree = degrees
    % Compute Future energy function
    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);

    fprintf('   %d   &   %8.2e   & %12.6e    \n', degree, tt, wzInit);
    futureTimes = [futureTimes, tt]; futureEnergies = [futureEnergies, wzInit];
end

% For run-time, only run the higher cases if exporting data, and only
% run once since the run time is longer so error is less sensitive
if exportData
    nTest = 1;
    for degree = 5:6
        % Compute Future energy function
        tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;

        wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);

        fprintf('   %d   &   %8.2e   & %12.6e    \n', degree, tt, wzInit);
        degrees = [degrees, degree]; futureTimes = [futureTimes, tt]; futureEnergies = [futureEnergies, wzInit];
    end
end
%% Export data
if exportData
    if x0 == 0.01
        fileName = sprintf('plots/example6_convergenceData_e%d.dat', numEl);
    else
        fileName = sprintf('plots/example6_convergenceData_e%d_biggerIC.dat', numEl);
    end
    fileID = fopen(fileName, 'w'); fprintf("Writing data to " + fileName + '\n')
    
    fprintf(fileID, '# Table III Data\n# finite element beam model, convergence and scalability results \n# numEls = %d   -->   n = %d \n   d   & CPU-sec-2    & E_d^+(x_0)      \n', numEl, 6 * numEl);
    
    for i = 1:length(degrees)
        fprintf(fileID, '   %d   &   %8.2e      & %12.6e    \n', degrees(i), futureTimes(i), futureEnergies(i));
    end
    fclose(fileID);
end


end
