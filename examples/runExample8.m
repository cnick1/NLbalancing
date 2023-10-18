function [w] = runExample8(exportData, x0)
%runExample8 Runs the finite element heat equation example to demonstrate
%            convergence and scalability.
%
%   Usage:  [v,w] = runExample8(plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportData)
%
%   Inputs:
%       exportData      - Boolean variable to determine if
%                         plots/data are exported
%       x0              - Initial condition
%
%   Outputs:
%       v,w             - Coefficients of the past and future energy
%                         function approximations, respectively
%
%   The value of eta is set below.
%
%   Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial drift nonlinear systems,” 2023.
%              [2] M. Embree, “Unstable modes in projection-based
%               reduced-order models: how many can there be, and what do they
%               tell you?,” Systems & Control Letters, vol. 124, pp. 49–59,
%               Feb. 2019, doi: 10.1016/j.sysconle.2018.11.010
%              [3] J. Galkowski, “Nonlinear instability in a semiclassical
%               problem,” Communications in Mathematical Physics, vol. 316,
%               no. 3, pp. 705–722, Oct. 2012, doi: 10.1007/s00220-012-1598-5
%              [4] B. Sandstede and A. Scheel, “Basin boundaries and
%               bifurcations near convective instabilities: a case study,”
%               Journal of Differential Equations, vol. 208, no. 1, pp.
%               176–193, Jan. 2005, doi: 10.1016/j.jde.2004.02.016
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

    [f, g, h] = getSystem8(numEl);

    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;

    fprintf(fileID, '%10.4e    & ', length(w{degree}));
    nd = [nd, length(w{degree})];
    fprintf(fileID, '%8.2e  & ', tt);
    times = [times, tt];

    % Initial condition from Mark Embree's talk
    L = 30; x = linspace(0, L, numEl + 1).';
    initialCondition = x0 * x .* (x - L) .* (x - L / 2);

    initialCondition = initialCondition(2:end - 1);

    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    fprintf(fileID, '%12.6e    \n', wzInit);
    energies = [energies, wzInit];
end

% For run-time, only run the higher cases if exporting data, and only
% run once since the run time is longer so error is less sensitive
if exportData
    nTest = 1;
    for numEl = [256, 512, 1024]
        numEls = [numEls, numEl];
        fprintf(fileID, '%5d       &', numEl); fprintf(fileID, '%5d & ', numEl - 1);

        [f, g, h] = getSystem8(numEl);

        tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;

        fprintf(fileID, '%10.4e    & ', length(w{degree}));
        nd = [nd, length(w{degree})];
        fprintf(fileID, '%8.2e  & ', tt);
        times = [times, tt];

        % Initial condition from Mark Embree's talk
        L = 30; x = linspace(0, L, numEl + 1).';
        initialCondition = x0 * x .* (x - L) .* (x - L / 2);

        initialCondition = initialCondition(2:end - 1);

        wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
        fprintf(fileID, '%12.6e    \n', wzInit);
        energies = [energies, wzInit];
    end
end
if exportData
    logy = log10(times); % take the natural log of y data
    logx = log10(numEls - 1); % take the natural log of x data
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
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e   &  %2.2f  &  %2.2f \n', numEls(i), numEls(i) - 1, nd(i), times(i), energies(i), a, d);
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

    [f, g, h] = getSystem8(numEl);

    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;

    fprintf(fileID, '%10.4e    & ', length(w{degree}));
    nd = [nd, length(w{degree})];
    fprintf(fileID, '%8.2e  & ', tt);
    times = [times, tt];

    % Initial condition from Mark Embree's talk
    L = 30; x = linspace(0, L, numEl + 1).';
    initialCondition = x0 * x .* (x - L) .* (x - L / 2);

    initialCondition = initialCondition(2:end - 1);

    wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
    fprintf(fileID, '%12.6e    \n', wzInit);
    energies = [energies, wzInit];
end

% For run-time, only run the higher cases if exporting data, and only
% run once since the run time is longer so error is less sensitive
if exportData
    nTest = 1;
    for numEl = [64, 128]
        numEls = [numEls, numEl];

        fprintf(fileID, '%5d       &', numEl); fprintf(fileID, '%5d & ', numEl - 1);

        [f, g, h] = getSystem8(numEl);

        tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;

        fprintf(fileID, '%10.4e    & ', length(w{degree}));
        nd = [nd, length(w{degree})];
        fprintf(fileID, '%8.2e  & ', tt);
        times = [times, tt];

        % Initial condition from Mark Embree's talk
        L = 30; x = linspace(0, L, numEl + 1).';
        initialCondition = x0 * x .* (x - L) .* (x - L / 2);

        initialCondition = initialCondition(2:end - 1);

        wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
        fprintf(fileID, '%12.6e    \n', wzInit);
        energies = [energies, wzInit];
    end
end

if exportData
    logy = log10(times); % take the natural log of y data
    logx = log10(numEls - 1); % take the natural log of x data
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
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e   &  %2.2f  &  %2.2f \n', numEls(i), numEls(i) - 1, nd(i), times(i), energies(i), a, d);
    end
    fclose(fileID);
end

% %%
% %  Computational performance of the energy function approximations.
% %  Since the initial times are so short, we average nTest times
% %
% %  This builds TABLE III
% %
% numEl = 32;
%
% fileID = 1; % Standard command window output if not writing to a file
%
% fprintf(fileID, '# Table III Data\n');
% fprintf(fileID, '# finite element heat equation model, convergence and scalability results \n');
% fprintf(fileID, '# numEls = %d   -->   n = %d \n', numEl, numEl-1);
%
% %print the header
% fprintf(fileID, 'd      ');
% % fprintf(fileID, '& CPU-sec   & E_d^-(x_0)     ');
% fprintf(fileID, '& CPU-sec-2    & E_d^+(x_0)      \n');
%
% % compute and print the results
% nTest = 3;
%
% [f, g, h] = getSystem8(numEl);
% % Initial condition where the nodes are displaced but have no initial
% % velocity or "rotation"
% numNodes = numEl + 1;
% % Initial condition from Mark Embree's talk
% L = 30; x = linspace(0,L,numEl+1).';
% initialCondition = x0 * x .* (x-L).*(x-L/2);
%
% initialCondition = initialCondition(2:end-1);
%
% pastTimes = []; futureTimes = []; pastEnergies = []; futureEnergies = [];
% degrees = 2:4;
% for degree = degrees
%     fprintf(fileID, '%d      & ', degree);
%
%     %     % Past
%     %     tic; for i = 1:nTest,
%     %         [v] = approxPastEnergy(f, g, C, eta, degree);
%     %     end, tt = toc / nTest;
%     %
%     %     fprintf(fileID, '%8.2e  & ', tt);
%     %     pastTimes = [pastTimes, tt];
%
%     %
%     %     vzInit = 0.5 * kronPolyEval(v, initialCondition, degree);
%     %     fprintf(fileID, '%12.6e    ', vzInit);
%     %     pastEnergies = [pastEnergies, vzInit];
%
%     % Future
%     tic; for i = 1:nTest,
%         [w] = approxFutureEnergy(f, g, h, eta, degree);
%     end, tt = toc / nTest;
%
%     fprintf(fileID, '%8.2e  & ', tt);
%     futureTimes = [futureTimes, tt];
%
%     wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
%     fprintf(fileID, '%12.6e    \n', wzInit);
%     futureEnergies = [futureEnergies, wzInit];
% end
%
% % For run-time, only run the higher cases if exporting data, and only
% % run once since the run time is longer so error is less sensitive
% if exportData
%     nTest = 1;
%     for degree = 5:6
%         degrees = [degrees, degree];
%         fprintf(fileID, '%d      & ', degree);
%
%         %     % Past
%         %     tic; for i = 1:nTest,
%         %         [v] = approxPastEnergy(f, g, C, eta, degree);
%         %     end, tt = toc / nTest;
%         %
%         %     fprintf(fileID, '%8.2e  & ', tt);
%         %     pastTimes = [pastTimes, tt];
%
%         %
%         %     vzInit = 0.5 * kronPolyEval(v, initialCondition, degree);
%         %     fprintf(fileID, '%12.6e    ', vzInit);
%         %     pastEnergies = [pastEnergies, vzInit];
%
%         % Future
%         tic; for i = 1:nTest,
%             [w] = approxFutureEnergy(f, g, h, eta, degree);
%         end, tt = toc / nTest;
%
%         fprintf(fileID, '%8.2e  & ', tt);
%         futureTimes = [futureTimes, tt];
%
%         wzInit = 0.5 * kronPolyEval(w, initialCondition, degree);
%         fprintf(fileID, '%12.6e    \n', wzInit);
%         futureEnergies = [futureEnergies, wzInit];
%
%     end
% end
% %% Export data
% if exportData
%     if x0 == 0.01
%         fileName = sprintf('plots/example8_convergenceData_e%d.dat', numEl);
%     else
%         fileName = sprintf('plots/example8_convergenceData_e%d_biggerIC.dat', numEl);
%     end
%     fprintf("Writing data to " + fileName + '\n')
%     fileID = fopen(fileName, 'w');
%
%     fprintf(fileID, '# Table III Data\n');
%     fprintf(fileID, '# finite element heat equation model, convergence and scalability results \n');
%     fprintf(fileID, '# numEls = %d   -->   n = %d \n', numEl, numEl-1);
%
%     %print the header
%     fprintf(fileID, 'd      ');
%     % fprintf(fileID, '& CPU-sec   & E_d^-(x_0)     ');
%     fprintf(fileID, '& CPU-sec-2    & E_d^+(x_0)      \n');
%     for i = 1:length(degrees)
%         fprintf(fileID, '%d      & ', degrees(i));
%         %     fprintf(fileID, '%8.2e  & ', pastTimes(i));
%         %     fprintf(fileID, '%12.6e    ', pastEnergies(i));
%         fprintf(fileID, '%8.2e     & ', futureTimes(i));
%         fprintf(fileID, '%12.6e    \n', futureEnergies(i));
%     end
%     fclose(fileID);
% end

end
