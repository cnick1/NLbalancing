function runExample6(exportData, x0)
%runExample6 Runs the finite element beam example to demonstrate
%            convergence and scalability.
%
%   Usage:  runExample6(exportData, x0)
%
%   Inputs:
%       exportData - Boolean to determine if data are exported
%       x0         - Initial condition scaling factor
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

eta = 0;

%% Table II: convergence & scalability wrt n (d=3)
degree = 3;

% print results to command window as we run, then write to file
fprintf('\n# Table II Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &   n  &     n^%d       &  CPU-sec  &  E_%d^+(x_0)     \n', degree, degree, degree);

nTest = [10, 10, 10, 10, 10, 1, 1, 1, 1];
numEls = [1, 2, 4, 8, 16, 32, 64, 128, 180];
nd = zeros(size(numEls)); times = zeros(size(numEls)); energies = zeros(size(numEls));
if exportData, numTestCases = 9; else, numTestCases = 5; end % If running with 16GB ram, 7; if 256GB ram, 9
for i=1:numTestCases
    fprintf(' %5d      &%4d  & ', numEls(i), 6 * numEls(i));
    
    % Compute energy functions and CPU time
    [f, g, h, IC] = getSystem6(numEls(i), 2);
    tic; for j = 1:nTest(i), [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest(i);
    
    % Evaluate energy function at x0 corresponding to nodes having linear
    % displacement but no rotation or initial velocity
    initialCondition = x0 * IC;
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    
    % Print results to command window
    fprintf(' %10.4e   & %8.2e  & %12.6e    \n', length(w{degree}), tt, wzInit);
    
    % log results to write to file later
    nd(i)= length(w{degree}); times(i) = tt; energies(i) = wzInit;
end

if exportData
    [a, d] = getExpFit(times, numEls);
    fprintf('The exponent fit gives n^%f compared with n^%d. \nWriting data to plots/example6_convergenceData_d3.dat \n', d, degree)
    fileID = fopen('plots/example6_convergenceData_d3.dat', 'w');
    fprintf(fileID, '# Table II Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &   n  &     n^%d       &  CPU-sec  &  E_%d^+(x_0)     &  exponentCoeff  &  exponentFit \n', degree, degree, degree);
    %print the header
    for i = 1:numTestCases
        fprintf(fileID, ' %5d      &%4d  &  %10.4e   & %8.2e  & %12.6e    &     %2.2f       &     %2.2f \n', numEls(i), 6 * numEls(i), nd(i), times(i), energies(i), a, d);
    end
    fclose(fileID);
end

%% Table II: convergence & scalability wrt n (d=4)
degree = 4;

%print the header
fprintf('\n# Table II Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &   n  &     n^%d       &  CPU-sec  &  E_%d^+(x_0)     \n', degree, degree, degree);

% compute and print the results
nTest = [3, 3, 3, 1, 1, 1];
numEls = [1, 2, 4, 8, 16, 24];
nd = zeros(size(numEls)); times = zeros(size(numEls)); energies = zeros(size(numEls)); 
if exportData, numTestCases = 6; else, numTestCases = 3; end
for i=1:numTestCases
    fprintf(' %5d      &%4d  & ', numEls(i), 6 * numEls(i));
    
    % Compute energy functions and CPU time
    [f, g, h, IC] = getSystem6(numEls(i), 2); 
    tic; for j = 1:nTest(i), [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest(i);
    
    % Evaluate energy function at x0 corresponding to nodes having linear
    % displacement but no rotation or initial velocity
    initialCondition = x0 * IC;
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    
    % Print results to command window
    fprintf(' %10.4e   & %8.2e  & %12.6e    \n', length(w{degree}), tt, wzInit);
    
    % log results to write to file later
    nd(i)= length(w{degree}); times(i) = tt; energies(i) = wzInit;
end
if exportData
    [a, d] = getExpFit(times, numEls);
    fprintf('The exponent fit gives n^%f compared with n^%d. \nWriting data to plots/example6_convergenceData_d4.dat \n', d, degree)
    fileID = fopen('plots/example6_convergenceData_d4.dat', 'w');
    fprintf(fileID, '# Table II Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &   n  &     n^%d       &  CPU-sec  &  E_%d^+(x_0)     &  exponentCoeff  &  exponentFit \n', degree, degree, degree);
    for i = 1:numTestCases
        fprintf(fileID, ' %5d      &%4d  &  %10.4e   & %8.2e  & %12.6e    &     %2.2f       &     %2.2f \n', numEls(i), 6 * numEls(i), nd(i), times(i), energies(i), a, d);
    end
    fclose(fileID);
end

%% Table III: convergence & scalability wrt d (n=18)
numEl = 3;

fprintf('\n# Table III Data\n# finite element beam model, convergence and scalability results \n# numEls = %d   -->   n = %d \n   d   &  CPU-sec-2   &  E_d^+(x_0)     \n', numEl, 6 * numEl);

% compute and print the results
[f, g, h, IC] = getSystem6(numEl, 2);

% Evaluate energy function at x0 corresponding to nodes having linear
% displacement but no rotation or initial velocity
initialCondition = x0 * IC;

degrees = 2:6;
nTest = [3, 3, 3, 3, 1];
futureTimes = zeros(size(degrees)); futureEnergies = zeros(size(degrees)); 
if exportData, numTestCases = 5; else, numTestCases = 3; end
for i=1:numTestCases
    % Compute Future energy function
    tic; for j = 1:nTest(i), [w] = approxFutureEnergy(f, g, h, eta, degrees(i)); end, tt = toc / nTest(i);
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degrees(i));

    fprintf('   %d   &   %8.2e   & %12.6e    \n', degrees(i), tt, wzInit);
    futureTimes(i) = tt; futureEnergies(i) = wzInit;
end

%% Export data
if exportData
    if x0 == 0.01
        fileName = sprintf('plots/example6_convergenceData_e%d.dat', numEl);
    else
        fileName = sprintf('plots/example6_convergenceData_e%d_biggerIC.dat', numEl);
    end
    fileID = fopen(fileName, 'w'); fprintf("Writing data to " + fileName + '\n')
    
    fprintf(fileID, '# Table III Data\n# finite element beam model, convergence and scalability results \n# numEls = %d   -->   n = %d \n   d   &  CPU-sec-2   &  E_d^+(x_0)     \n', numEl, 6 * numEl);
    
    for i = 1:numTestCases
        fprintf(fileID, '   %d   &   %8.2e   & %12.6e    \n', degrees(i), futureTimes(i), futureEnergies(i));
    end
    fclose(fileID);
end


end

function [a, d] = getExpFit(times, numEls)
n = find(times~=0,1,'last'); times = times(1:n); numEls = numEls(1:n); 
% CPU time scales as O(n^d) for some d; this return d
logy = log10(times); % take the natural log of y data
logx = log10(6 * numEls); % take the natural log of x data
X = [ones(length(logx), 1) logx']; % create matrix of x data with a column of ones
beta = X \ logy'; % solve for beta coefficients using linear least squares
a = beta(1); % line y intercept (moves line up and down)
d = beta(2); % calculate exponent d (slope of the line)
end
