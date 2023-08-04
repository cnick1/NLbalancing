function [w] = runExample8(exportData, x0, varargin)
%runExample8 Runs the finite element heat equation example to demonstrate
%            convergence and scalability.
%
%   Usage:  [v,w] = runExample8(plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportData)
%
%   Inputs:
%       exportData      - Boolean variable to determine if
%                         plots/data are exported
%
%   Outputs:
%       v,w             - Coefficients of the past and future energy
%                         function approximations, respectively
%
%   The value of eta is set below.
%
%   Part of the NLbalancing repository.
%%

if nargin < 2
    if nargin < 1
        exportData = false;
    end
    x0 = 1e-5;
end

eta = 0.5;

fprintf('Running Example 8\n')

%%
%  Computational performance of the energy function approximations.
%  Since the initial times are so short, we average nTest times
%
%  This builds TABLE I
%
fileID = 1;
degree = 3;

fprintf(fileID, '# Table I Data\n');
fprintf(fileID, '# finite element heat equation model, convergence and scalability results; d=%d \n', degree);
%print the header
fprintf(fileID, 'numElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)      \n', degree, degree);

% compute and print the results
nTest = 10;
nd = []; times = []; energies = [];
numEls = [4, 8, 16, 32, 64, 128];
for numEl = numEls
    fprintf(fileID, '%5d       &', numEl); fprintf(fileID, '%5d & ', numEl - 1);
    
    [A, B, C, N, f, g, h] = getSystem8(numEl);
    
    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, N, g, C, eta, degree); end, tt = toc / nTest;
    
    fprintf(fileID, '%10.4e    & ', length(w{degree}));
    nd = [nd, length(w{degree})];
    fprintf(fileID, '%8.2e  & ', tt);
    times = [times, tt];
    
    % Initial condition from Mark Embree's talk
    L = 30; x = linspace(0,L,numEl+1).';
    initialCondition = x0 * x .* (x-L).*(x-L/2);
    
    initialCondition = initialCondition(2:end-1);
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    fprintf(fileID, '%12.6e    \n', wzInit);
    energies = [energies, wzInit];
end

% For run-time, only run the higher cases if exporting data, and only
% run once since the run time is longer so error is less sensitive
if exportData
    nTest = 1;
    numEls = [256, 512, 1024, 1024];
    for numEl = numEls
        fprintf(fileID, '%5d       &', numEl); fprintf(fileID, '%5d & ', numEl - 1);
        
        [A, B, C, N, f, g, h] = getSystem8(numEl);
        
        tic; for i = 1:nTest, [w] = approxFutureEnergy(f, N, g, C, eta, degree); end, tt = toc / nTest;
        
        fprintf(fileID, '%10.4e    & ', length(w{degree}));
        nd = [nd, length(w{degree})];
        fprintf(fileID, '%8.2e  & ', tt);
        times = [times, tt];
        
        % Initial condition from Mark Embree's talk
        L = 30; x = linspace(0,L,numEl+1).';
        initialCondition = x0 * x .* (x-L).*(x-L/2);
                
        initialCondition = initialCondition(2:end-1);
        
        wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
        fprintf(fileID, '%12.6e    \n', wzInit);
        energies = [energies, wzInit];
    end
end
if exportData
    logy = log10(times); % take the natural log of y data
    logx = log10(numEls-1); % take the natural log of x data
    X = [ones(length(logx), 1) logx']; % create matrix of x data with a column of ones
    beta = X \ logy'; % solve for beta coefficients using linear least squares
    a = beta(1); % calculate exponent d
    d = beta(2); % calculate exponent d
    fprintf('The exponent fit gives n^%f compared with n^%d. \n', d, degree)
    fprintf('Writing data to plots/example8_convergenceData_d3.dat \n')
    fileID = fopen('plots/example8_convergenceData_d3.dat', 'w');
    fprintf(fileID, '# Table I Data\n');
    fprintf(fileID, '# finite element heat equation model, convergence and scalability results; d=%d \n', degree);
    %print the header
    fprintf(fileID, 'numElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)     &  exponentCoeff  &  exponentFit \n', degree, degree);
    for i = 1:length(numEls)
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e   &  %2.2f  &  %2.2f \n', numEls(i), numEls(i)-1, nd(i), times(i), energies(i), a, d);
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
fprintf(fileID, '# finite element heat equation model, convergence and scalability results; d=%d \n', degree);
%print the header
fprintf(fileID, 'numElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)      \n', degree, degree);

% compute and print the results
nTest = 3;
nd = []; times = []; energies = [];
numEls = [4, 8, 16, 32];
for numEl = numEls
    fprintf(fileID, '%5d       &', numEl); fprintf(fileID, '%5d & ', numEl - 1);
    
    [A, B, C, N, f, g, h] = getSystem8(numEl);
    
    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, N, g, C, eta, degree); end, tt = toc / nTest;
    
    fprintf(fileID, '%10.4e    & ', length(w{degree}));
    nd = [nd, length(w{degree})];
    fprintf(fileID, '%8.2e  & ', tt);
    times = [times, tt];
    
    % Initial condition from Mark Embree's talk
    L = 30; x = linspace(0,L,numEl+1).';
    initialCondition = x0 * x .* (x-L).*(x-L/2);
        
    initialCondition = initialCondition(2:end-1);
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    fprintf(fileID, '%12.6e    \n', wzInit);
    energies = [energies, wzInit];
end

% For run-time, only run the higher cases if exporting data, and only
% run once since the run time is longer so error is less sensitive
if exportData
    nTest = 1;
    numEls = [64, 128];
    for numEl = numEls
        fprintf(fileID, '%5d       &', numEl); fprintf(fileID, '%5d & ', numEl - 1);
        
        [A, B, C, N, f, g, h] = getSystem8(numEl);
        
        tic; for i = 1:nTest, [w] = approxFutureEnergy(f, N, g, C, eta, degree); end, tt = toc / nTest;
        
        fprintf(fileID, '%10.4e    & ', length(w{degree}));
        nd = [nd, length(w{degree})];
        fprintf(fileID, '%8.2e  & ', tt);
        times = [times, tt];
        
        % Initial condition from Mark Embree's talk
        L = 30; x = linspace(0,L,numEl+1).';
        initialCondition = x0 * x .* (x-L).*(x-L/2);
        
        initialCondition = initialCondition(2:end-1);
        
        wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
        fprintf(fileID, '%12.6e    \n', wzInit);
        energies = [energies, wzInit];
    end
end

if exportData
    logy = log10(times); % take the natural log of y data
    logx = log10(numEls-1); % take the natural log of x data
    X = [ones(length(logx), 1) logx']; % create matrix of x data with a column of ones
    beta = X \ logy'; % solve for beta coefficients using linear least squares
    a = beta(1); % calculate exponent d
    d = beta(2); % calculate exponent d
    fprintf('The exponent fit gives n^%f compared with n^%d. \n', d, degree)
    fprintf('Writing data to plots/example8_convergenceData_d4.dat \n')
    fileID = fopen('plots/example8_convergenceData_d4.dat', 'w');
    fprintf(fileID, '# Table II Data\n');
    fprintf(fileID, '# finite element heat equation model, convergence and scalability results; d=%d \n', degree);
    %print the header
    fprintf(fileID, 'numElements &    n & n^%d           & CPU-sec   & E_%d^+(x_0)     &  exponentCoeff  &  exponentFit \n', degree, degree);
    for i = 1:length(numEls)
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e   &  %2.2f  &  %2.2f \n', numEls(i), numEls(i)-1, nd(i), times(i), energies(i), a, d);
    end
    fclose(fileID);
end

%%
%  Computational performance of the energy function approximations.
%  Since the initial times are so short, we average nTest times
%
%  This builds TABLE III
%
numEl = 32;

fileID = 1; % Standard command window output if not writing to a file

fprintf(fileID, '# Table III Data\n');
fprintf(fileID, '# finite element heat equation model, convergence and scalability results \n');
fprintf(fileID, '# numEls = %d   -->   n = %d \n', numEl, numEl-1);

%print the header
fprintf(fileID, 'd      ');
% fprintf(fileID, '& CPU-sec   & E_d^-(x_0)     ');
fprintf(fileID, '& CPU-sec-2    & E_d^+(x_0)      \n');

% compute and print the results
nTest = 3;

[A, B, C, N, f, g, h] = getSystem8(numEl);
% Initial condition where the nodes are displaced but have no initial
% velocity or "rotation"
numNodes = numEl + 1;
% Initial condition from Mark Embree's talk
L = 30; x = linspace(0,L,numEl+1).';
initialCondition = x0 * x .* (x-L).*(x-L/2);

initialCondition = initialCondition(2:end-1);

pastTimes = []; futureTimes = []; pastEnergies = []; futureEnergies = [];
degrees = 2:4;
for degree = degrees
    fprintf(fileID, '%d      & ', degree);
    
    %     % Past
    %     tic; for i = 1:nTest,
    %         [v] = approxPastEnergy(A, N, g, C, eta, degree);
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
        [w] = approxFutureEnergy(f, N, g, C, eta, degree);
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
        %         [v] = approxPastEnergy(A, N, g, C, eta, degree);
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
            [w] = approxFutureEnergy(f, N, g, C, eta, degree);
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
        fileName = sprintf('plots/example8_convergenceData_e%d.dat', numEl);
    else
        fileName = sprintf('plots/example8_convergenceData_e%d_biggerIC.dat', numEl);
    end
    fprintf("Writing data to " + fileName + '\n')
    fileID = fopen(fileName, 'w');
    
    fprintf(fileID, '# Table III Data\n');
    fprintf(fileID, '# finite element heat equation model, convergence and scalability results \n');
    fprintf(fileID, '# numEls = %d   -->   n = %d \n', numEl, numEl-1);
    
    %print the header
    fprintf(fileID, 'd      ');
    % fprintf(fileID, '& CPU-sec   & E_d^-(x_0)     ');
    fprintf(fileID, '& CPU-sec-2    & E_d^+(x_0)      \n');
    for i = 1:length(degrees)
        fprintf(fileID, '%d      & ', degrees(i));
        %     fprintf(fileID, '%8.2e  & ', pastTimes(i));
        %     fprintf(fileID, '%12.6e    ', pastEnergies(i));
        fprintf(fileID, '%8.2e     & ', futureTimes(i));
        fprintf(fileID, '%12.6e    \n', futureEnergies(i));
    end
    fclose(fileID);
end


end
