function [w] = runExample6(numGTermsModel, numGTermsApprox, exportData, x0, varargin)
%runExample6 Runs the finite element beam example to demonstrate
%            convergence and scalability.
%
%   Usage:  [v,w] = runExample6(plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportData)
%
%   Inputs:
%       numGTermsModel  - Number of terms in the full order model
%       numGTermsApprox - Number of terms assumed when computing energy
%                         functions
%       exportData      - Boolean variable to determine if
%                         plots/data are exported
%
%   Outputs:
%       v,w             - Coefficients of the past and future energy
%                         function approximations, respectively
%
%   The value of eta is set below.
%
%   Reference: Submitted to IEEE TAC. 
%   "Scalable Computation of ℋ∞ Energy Functions for Polynomial 
%    Control-Affine Systems", N. Corbin and B. Kramer
%    arXiv:
%
%   Part of the NLbalancing repository.
%%

if nargin < 4
    if nargin < 3
        if nargin < 2
            if nargin < 1
                numGTermsModel = 2;
            end
            numGTermsApprox = numGTermsModel;
        end
        exportData = false;
    end
    x0 = 0.01;
end

eta = 0;

fprintf('Running Example 6\n')

%%
%  Computational performance of the energy function approximations.
%  Since the initial times are so short, we average nTest times
%
%  This builds TABLE I
%
fileID = 1;
degree = 3;


fprintf(fileID, '# Table I Data\n');
fprintf(fileID, '# finite element beam model, convergence and scalability results; d=%d \n', degree);
%print the header
fprintf(fileID, 'numElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)      \n', degree, degree);

% compute and print the results
nTest = 10;
nd = []; times = []; energies = [];
numEls = [1, 2, 4, 8, 16];
for numEl = numEls
    fprintf(fileID, '%5d       &', numEl);
    fprintf(fileID, '%5d & ', 6 * numEl);
    [A, B, C, N, f, g, h] = getSystem6(numEl,2);
    f{4} = sparse(length(A),length(A)^4);
    f = f(1:numGTermsModel); % Adjust FOM to be Quadratic, QB, etc.
    g = g(1:numGTermsModel); % Adjust FOM to be Quadratic, QB, etc.
    
    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, N, g(1:numGTermsApprox), C, eta, degree); end, tt = toc / nTest;
    
    fprintf(fileID, '%10.4e    & ', length(w{degree}));
    nd = [nd, length(w{degree})];
    fprintf(fileID, '%8.2e  & ', tt);
    times = [times, tt];
    
    % Initial condition where the nodes are displaced but have no initial
    % velocity or "rotation"
    numNodes = numEl + 1;
    initialCondition = x0 / (numNodes - 1) * ...
        [[(0:numNodes - 1);
        (0:numNodes - 1);
        0 * (0:numNodes - 1)].';
        zeros(numNodes, 3)].'; % Full initial condition
    initialCondition(:, 1) = []; initialCondition(:, 1 + numNodes) = []; % Remove the first node DOFs
    initialCondition = initialCondition(:);
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    fprintf(fileID, '%12.6e    \n', wzInit);
    energies = [energies, wzInit];
end

% For run-time, only run the higher cases if exporting data, and only
% run once since the run time is longer so error is less sensitive
if exportData
    nTest = 1;
    for numEl = [32, 64, 128] %128 runs out of ram in kroneckerLeft.m
        numEls = [numEls, numEl];
        fprintf(fileID, '%5d       &', numEl);
        fprintf(fileID, '%5d & ', 6 * numEl);
        [A, B, C, N, f, g, h] = getSystem6(numEl);
%     f{4} = sparse(length(A),length(A)^4);
%     f = f(1:numGTermsModel); % Adjust FOM to be Quadratic, QB, etc.
    g = g(1:numGTermsModel); % Adjust FOM to be Quadratic, QB, etc.        
        tic; for i = 1:nTest, [w] = approxFutureEnergy(f, N, g(1:numGTermsApprox), C, eta, degree); end, tt = toc / nTest;
        
        fprintf(fileID, '%10.4e    & ', length(w{degree}));
        nd = [nd, length(w{degree})];
        fprintf(fileID, '%8.2e  & ', tt);
        times = [times, tt];
        
        % Initial condition where the nodes are displaced but have no initial
        % velocity or "rotation"
        numNodes = numEl + 1;
        initialCondition = x0 / (numNodes - 1) * ...
            [[(0:numNodes - 1);
            (0:numNodes - 1);
            0 * (0:numNodes - 1)].';
            zeros(numNodes, 3)].'; % Full initial condition
        initialCondition(:, 1) = []; initialCondition(:, 1 + numNodes) = []; % Remove the first node DOFs
        initialCondition = initialCondition(:);
        
        wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
        fprintf(fileID, '%12.6e    \n', wzInit);
        energies = [energies, wzInit];
    end
end
if false%exportData
    logy = log10(times); % take the natural log of y data
    logx = log10(6*numEls); % take the natural log of x data
    X = [ones(length(logx),1) logx']; % create matrix of x data with a column of ones
    beta = X\logy'; % solve for beta coefficients using linear least squares
    a = beta(1); % calculate exponent d
    d = beta(2); % calculate exponent d
    fprintf('The exponent fit gives n^%f compared with n^%d. \n', d, degree)
    fprintf('Writing data to plots/example6_convergenceData_d3.dat \n')
    fileID = fopen('plots/example6_convergenceData_d3.dat', 'w');
    fprintf(fileID, '# Table I Data\n');
    fprintf(fileID, '# finite element beam model, convergence and scalability results; d=%d \n', degree);
    %print the header
    fprintf(fileID, 'numElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)     &  exponentCoeff  &  exponentFit \n', degree, degree);
    for i=1:length(numEls)
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e   &  %2.2f  &  %2.2f \n', numEls(i),6*numEls(i),nd(i),times(i),energies(i), a, d);
    end
    fclose(fileID);
end

%%
%  Computational performance of the energy function approximations.
%  Since the initial times are so short, we average nTest times
%
%  This builds TABLE II
%
fileID = 1;
degree = 4;


fprintf(fileID, '# Table II Data\n');
fprintf(fileID, '# finite element beam model, convergence and scalability results; d=%d \n', degree);
%print the header
fprintf(fileID, 'numElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)      \n', degree, degree);

% compute and print the results
nTest = 3;
nd = []; times = []; energies = [];
numEls = [1, 2, 4];
for numEl = numEls
    fprintf(fileID, '%5d       &', numEl);
    fprintf(fileID, '%5d & ', 6 * numEl);
    [A, B, C, N, f, g, h] = getSystem6(numEl);
%     f{4} = sparse(length(A),length(A)^4);
%     f = f(1:numGTermsModel); % Adjust FOM to be Quadratic, QB, etc.
    g = g(1:numGTermsModel); % Adjust FOM to be Quadratic, QB, etc.    
    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, N, g(1:numGTermsApprox), C, eta, degree); end, tt = toc / nTest;
    
    fprintf(fileID, '%10.4e    & ', length(w{degree}));
    nd = [nd, length(w{degree})];
    fprintf(fileID, '%8.2e  & ', tt);
    times = [times, tt];
    
    % Initial condition where the nodes are displaced but have no initial
    % velocity or "rotation"
    numNodes = numEl + 1;
    initialCondition = x0 / (numNodes - 1) * ...
        [[(0:numNodes - 1);
        (0:numNodes - 1);
        0 * (0:numNodes - 1)].';
        zeros(numNodes, 3)].'; % Full initial condition
    initialCondition(:, 1) = []; initialCondition(:, 1 + numNodes) = []; % Remove the first node DOFs
    initialCondition = initialCondition(:);
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    fprintf(fileID, '%12.6e    \n', wzInit);
    energies = [energies, wzInit];
end

% For run-time, only run the higher cases if exporting data, and only
% run once since the run time is longer so error is less sensitive
if exportData
    nTest = 1;
    for numEl = [8, 16, 32, 64]
        numEls = [numEls, numEl];
        fprintf(fileID, '%5d       &', numEl);
        fprintf(fileID, '%5d & ', 6 * numEl);
        [A, B, C, N, f, g, h] = getSystem6(numEl);
    f{4} = sparse(length(A),length(A)^4);
    f = f(1:numGTermsModel); % Adjust FOM to be Quadratic, QB, etc.
    g = g(1:numGTermsModel); % Adjust FOM to be Quadratic, QB, etc.        
        tic; for i = 1:nTest, [w] = approxFutureEnergy(f, N, g(1:numGTermsApprox), C, eta, degree); end, tt = toc / nTest;
        
        fprintf(fileID, '%10.4e    & ', length(w{degree}));
        nd = [nd, length(w{degree})];
        fprintf(fileID, '%8.2e  & ', tt);
        times = [times, tt];
        
        % Initial condition where the nodes are displaced but have no initial
        % velocity or "rotation"
        numNodes = numEl + 1;
        initialCondition = x0 / (numNodes - 1) * ...
            [[(0:numNodes - 1);
            (0:numNodes - 1);
            0 * (0:numNodes - 1)].';
            zeros(numNodes, 3)].'; % Full initial condition
        initialCondition(:, 1) = []; initialCondition(:, 1 + numNodes) = []; % Remove the first node DOFs
        initialCondition = initialCondition(:);
        
        wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
        fprintf(fileID, '%12.6e    \n', wzInit);
        energies = [energies, wzInit];
    end
end
if exportData
    logy = log10(times); % take the natural log of y data
    logx = log10(6*numEls); % take the natural log of x data
    X = [ones(length(logx),1) logx']; % create matrix of x data with a column of ones
    beta = X\logy'; % solve for beta coefficients using linear least squares
    a = beta(1); % calculate exponent d
    d = beta(2); % calculate exponent d
    fprintf('The exponent fit gives n^%f compared with n^%d. \n', d, degree)
    fprintf('Writing data to plots/example6_convergenceData_d4.dat \n')
    fileID = fopen('plots/example6_convergenceData_d4.dat', 'w');
    fprintf(fileID, '# Table II Data\n');
    fprintf(fileID, '# finite element beam model, convergence and scalability results; d=%d \n', degree);
    %print the header
    fprintf(fileID, 'numElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)     &  exponentCoeff  &  exponentFit \n', degree, degree);
    for i=1:length(numEls)
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e   &  %2.2f  &  %2.2f \n', numEls(i),6*numEls(i),nd(i),times(i),energies(i), a, d);
    end
    fclose(fileID);
end


%%
%  Computational performance of the energy function approximations.
%  Since the initial times are so short, we average nTest times
%
%  This builds TABLE III
%
numEl = 3;

fileID = 1; % Standard command window output if not writing to a file


fprintf(fileID, '# Table III Data\n');
fprintf(fileID, '# finite element beam model, convergence and scalability results \n');
fprintf(fileID, '# numEls = %d   -->   n = %d \n', numEl, 6*numEl);

%print the header
fprintf(fileID, 'd      ');
% fprintf(fileID, '& CPU-sec   & E_d^-(x_0)     ');
fprintf(fileID, '& CPU-sec-2    & E_d^+(x_0)      \n');

% compute and print the results
nTest = 3;

[A, B, C, N, f, g, h] = getSystem6(numEl);
    f{4} = sparse(length(A),length(A)^4);
    f = f(1:numGTermsModel); % Adjust FOM to be Quadratic, QB, etc.
    g = g(1:numGTermsModel); % Adjust FOM to be Quadratic, QB, etc.
% Initial condition where the nodes are displaced but have no initial
% velocity or "rotation"
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
    fprintf(fileID, '%d      & ', degree);
    
    %     % Past
    %     tic; for i = 1:nTest,
    %         [v] = approxPastEnergy(A, N, g(1:numGTermsApprox), C, eta, degree);
    %     end, tt = toc / nTest;
    %
    %     fprintf(fileID, '%8.2e  & ', tt);
    %     pastTimes = [pastTimes, tt];
    
    %
    %     vzInit = 0.5 * kronPolyEval(v, initialCondition, degree);
    %     fprintf(fileID, '%12.6e    ', vzInit);
    %     pastEnergies = [pastEnergies, vzInit];
    
    
    % Future
    tic; for i = 1:nTest,
        [w] = approxFutureEnergy(f, N, g(1:numGTermsApprox), C, eta, degree);
    end, tt = toc / nTest;
    
    fprintf(fileID, '%8.2e  & ', tt);
    futureTimes = [futureTimes, tt];
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    fprintf(fileID, '%12.6e    \n', wzInit);
    futureEnergies = [futureEnergies, wzInit];
end

% For run-time, only run the higher cases if exporting data, and only
% run once since the run time is longer so error is less sensitive
if exportData
    nTest = 1;
    for degree = 5:6
        degrees = [degrees, degree];
        fprintf(fileID, '%d      & ', degree);
        
        %     % Past
        %     tic; for i = 1:nTest,
        %         [v] = approxPastEnergy(A, N, g(1:numGTermsApprox), C, eta, degree);
        %     end, tt = toc / nTest;
        %
        %     fprintf(fileID, '%8.2e  & ', tt);
        %     pastTimes = [pastTimes, tt];
        
        %
        %     vzInit = 0.5 * kronPolyEval(v, initialCondition, degree);
        %     fprintf(fileID, '%12.6e    ', vzInit);
        %     pastEnergies = [pastEnergies, vzInit];
        
        
        % Future
        tic; for i = 1:nTest,
            [w] = approxFutureEnergy(f, N, g(1:numGTermsApprox), C, eta, degree);
        end, tt = toc / nTest;
        
        fprintf(fileID, '%8.2e  & ', tt);
        futureTimes = [futureTimes, tt];
        
        wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
        fprintf(fileID, '%12.6e    \n', wzInit);
        futureEnergies = [futureEnergies, wzInit];
        
    end
end
%% Export data
if exportData
    if x0 == 0.01
    	fileName = sprintf('plots/example6_convergenceData_e%d.dat',numEl);
    else
    	fileName = sprintf('plots/example6_convergenceData_e%d_biggerIC.dat',numEl);
    end
    fprintf("Writing data to " + fileName + '\n')
    fileID = fopen(fileName, 'w');
    
    fprintf(fileID, '# Table III Data\n');
    fprintf(fileID, '# finite element beam model, convergence and scalability results \n');
    fprintf(fileID, '# numEls = %d   -->   n = %d \n', numEl, 6*numEl);
    
    %print the header
    fprintf(fileID, 'd      ');
    % fprintf(fileID, '& CPU-sec   & E_d^-(x_0)     ');
    fprintf(fileID, '& CPU-sec-2    & E_d^+(x_0)      \n');
    for i=1:length(degrees)
        fprintf(fileID, '%d      & ', degrees(i));
        %     fprintf(fileID, '%8.2e  & ', pastTimes(i));
        %     fprintf(fileID, '%12.6e    ', pastEnergies(i));
        fprintf(fileID, '%8.2e     & ', futureTimes(i));
        fprintf(fileID, '%12.6e    \n', futureEnergies(i));
    end
    fclose(fileID);
end
%% Perform curve fitting for finding n^alpha scaling
%
%
% opts = delimitedTextImportOptions("NumVariables", 5); opts.DataLines = [4, Inf]; opts.Delimiter = "&"; opts.VariableNames = ["numElements", "n", "n4", "CPUsec", "E_4x_0"]; opts.VariableTypes = ["double", "double", "double", "double", "double"]; opts.ExtraColumnsRule = "ignore"; opts.EmptyLineRule = "read";
% example6convergenceDatad4 = readtable("D:\share\NLbalancing\plots\example6_convergenceData_d4.dat", opts);
% clear opts
%
%
% logy = log10example6convergenceDatad4.CPUsec).'; % take the natural log of y data
% logx = log10example6convergenceDatad4.n).'; % take the natural log of x data
% X = [ones(length(logx),1) logx']; % create matrix of x data with a column of ones
% beta = X\logy'; % solve for beta coefficients using linear least squares
% d = beta(2) % calculate exponent d


end
