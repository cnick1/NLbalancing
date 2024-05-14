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

eta = 0.1;

%% Table II: convergence & scalability wrt n (d=3)
degree = 3;

% print results to command window as we run, then write to file
fprintf('\n# Table II Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &   n  &     n^%d       &  CPU-sec  &  E_%d^+(x_0)     &  RES_%d^+(x_0)     \n', degree, degree, degree, degree);

nTest = [10, 10, 10, 10, 10, 1, 1, 1, 1];
numEls = [1, 2, 4, 8, 16, 32, 64, 128, 180];
nd = zeros(size(numEls)); times = zeros(size(numEls)); energies = zeros(size(numEls)); RES = zeros(size(numEls));
if exportData, numTestCases = 6; else, numTestCases = 5; end % If running with 16GB ram, 7; if 256GB ram, 9
for i=1:numTestCases
    fprintf(' %5d      &%4d  & ', numEls(i), 6 * numEls(i));
    
    % Compute energy functions and CPU time
    [f, g, h, IC] = getSystem6(numEls(i), 2);
    tic; for j = 1:nTest(i), [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest(i);
    
    % Evaluate energy function at x0 corresponding to nodes having linear
    % displacement but no rotation or initial velocity
    initialCondition = x0 * IC;
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    
    % Evaluate HJB residual
    RES(i) = getFutureHJBResidual(f, g, h, eta, w, initialCondition);
    
    % Print results to command window
    fprintf(' %10.4e   & %8.2e  & %12.6e    &  %12.6e      \n', length(w{degree}), tt, wzInit, RES(i));
    
    % log results to write to file later
    nd(i)= length(w{degree}); times(i) = tt; energies(i) = wzInit;
end

if exportData
    [a, d] = getExpFit(times, numEls);
    fprintf('The exponent fit gives n^%f compared with n^%d. \nWriting data to plots/example6_convergenceData_d3.dat \n', d, degree)
    fileID = fopen('plots/example6_convergenceData_d3.dat', 'w');
    fprintf(fileID, '# Table I Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)     &  RES_%d^+(x_0)     &  exponentCoeff  &  exponentFit \n', degree, degree, degree, degree);
    %print the header
    for i = 1:numTestCases
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e    &  %12.6e      &  %2.2f  &  %2.2f \n', numEls(i), 6 * numEls(i), nd(i), times(i), energies(i), RES(i), a, d);
    end
    fclose(fileID);
end

%% Table II: convergence & scalability wrt n (d=4)
degree = 4;

%print the header
fprintf('\n# Table II Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &   n  &     n^%d       &  CPU-sec  &  E_%d^+(x_0)     &  RES_%d^+(x_0)     \n', degree, degree, degree, degree);

% compute and print the results
nTest = [3, 3, 3, 1, 1];
numEls = [1, 2, 4, 8, 16];
nd = zeros(size(numEls)); times = zeros(size(numEls)); energies = zeros(size(numEls)); RES = zeros(size(numEls));
if exportData, numTestCases = 4; else, numTestCases = 3; end
for i=1:numTestCases
    fprintf(' %5d      &%4d  & ', numEls(i), 6 * numEls(i));
    
    % Compute energy functions and CPU time
    [f, g, h, IC] = getSystem6(numEls(i), 2);
    tic; for j = 1:nTest(i), [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest(i);
    
    % Evaluate energy function at x0 corresponding to nodes having linear
    % displacement but no rotation or initial velocity
    initialCondition = x0 * IC;
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    
    % Evaluate HJB residual
    RES(i) = getFutureHJBResidual(f, g, h, eta, w, initialCondition);
    
    % Print results to command window
    fprintf(' %10.4e   & %8.2e  & %12.6e    &  %12.6e      \n', length(w{degree}), tt, wzInit, RES(i));
    
    % log results to write to file later
    nd(i)= length(w{degree}); times(i) = tt; energies(i) = wzInit;
end
if exportData
    [a, d] = getExpFit(times, numEls);
    fprintf('The exponent fit gives n^%f compared with n^%d. \nWriting data to plots/example6_convergenceData_d4.dat \n', d, degree)
    fileID = fopen('plots/example6_convergenceData_d4.dat', 'w');
    fprintf(fileID, '# Table II Data\n# finite element beam model, convergence and scalability results; d=%d \nnumElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)     &  RES_%d^+(x_0)     &  exponentCoeff  &  exponentFit \n', degree, degree, degree, degree);
    for i = 1:numTestCases
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e    &  %12.6e      &  %2.2f  &  %2.2f \n', numEls(i), 6 * numEls(i), nd(i), times(i), energies(i), RES, a, d);
    end
    fclose(fileID);
end

%% Table III: convergence & scalability wrt d (n=18)
numEl = 4;

fprintf('\n# Table III Data\n# finite element beam model, convergence and scalability results \n# numEls = %d   -->   n = %d \n   d   &  CPU-sec-2   & E_d^+(x_0)     &  RES_d^+(x_0)     \n', numEl, 6 * numEl);

% compute and print the results
[f, g, h, IC] = getSystem6(numEl, 2);

% Evaluate energy function at x0 corresponding to nodes having linear
% displacement but no rotation or initial velocity
initialCondition = x0 * IC;

degrees = 2:6;
nTest = [3, 3, 3, 3, 1];
futureTimes = zeros(size(degrees)); futureEnergies = zeros(size(degrees)); RES = zeros(size(degrees));
if exportData, numTestCases = 5; else, numTestCases = 3; end
for i=1:numTestCases
    % Compute Future energy function
    tic; for j = 1:nTest(i), [w] = approxFutureEnergy(f, g, h, eta, degrees(i)); end, tt = toc / nTest(i);
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degrees(i));
        
    % Evaluate HJB residual
    RES(i) = getFutureHJBResidual(f, g, h, eta, w, initialCondition);

    fprintf('   %d   &   %8.2e   & %12.6e    & %12.6e    \n', degrees(i), tt, wzInit, RES(i));
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
    
    fprintf(fileID, '# Table III Data\n# finite element beam model, convergence and scalability results \n# numEls = %d   -->   n = %d \n   d   & CPU-sec-2    &  E_d^+(x_0)     \n', numEl, 6 * numEl);
    
    for i = 1:numTestCases
        fprintf(fileID, '   %d   &   %8.2e      & %12.6e    & %12.6e    \n', degrees(i), futureTimes(i), futureEnergies(i), RES(i));
    end
    fclose(fileID);
end


end

function [a, d] = getExpFit(times, numEls)
% CPU time scales as O(n^d) for some d; this return d
logy = log10(times); % take the natural log of y data
logx = log10(6 * numEls); % take the natural log of x data
X = [ones(length(logx), 1) logx']; % create matrix of x data with a column of ones
beta = X \ logy'; % solve for beta coefficients using linear least squares
a = beta(1); % line y intercept (moves line up and down)
d = beta(2); % calculate exponent d (slope of the line)
end

function RES = getFutureHJBResidual(f, g, h, eta, w, x)
% Evaluate the Future HJB PDE at the point x; should be zero, so whatever
% it is represents the error (residual)
dEdX = 0.5 * kronPolyDerivEval(w, x);
FofX = kronPolyEval(f, x);
GofX = g{1}; lg = length(g); Im = speye(size(g{1},2));
xk = 1; for k=2:lg; xk = kron(xk, x); GofX = GofX + g{k} * kron(xk, Im); end
HofX = kronPolyEval(h, x);

RES = dEdX*FofX - eta/2*dEdX*(GofX*GofX.')*dEdX.' + 0.5*(HofX.'*HofX);
end