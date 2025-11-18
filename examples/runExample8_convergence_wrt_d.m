function runExample8_convergence_wrt_d(exportData, x0)
%runExample8_convergence_wrt_d Runs the finite element heat equation example to demonstrate
%            convergence and scalability.
%
%   Usage:  runExample8_convergence_wrt_d(plotEnergy,plotBalancing,balancingDegree,
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
%   Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of
%               𝓗∞ energy functions for polynomial drift nonlinear
%               systems,” in 2024 American Control Conference (ACC), Jul.
%               2024, pp. 2506–2511. doi: 10.23919/acc60939.2024.10644363
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
    x0 = 5e-5;
end

eta = 0.5;

fprintf('Running Example 8\n')

%%
%  Computational performance of the energy function approximations.
%  This builds TABLE II
%

%print the header
fprintf('# Table II Data\n# finite element heat equation model, convergence and scalability results; \nd & E_d^+(x_0)      \n');

% compute and print the results
energies = [];
if exportData
    numEl = 16;
    ds = 2:2:6;
else
    numEl = 8;
    ds = 2:2:8;
end
for d = ds
    fprintf('%1d & ', d);
    
    [f, g, h] = getSystem8(numEl);
    
    [w] = approxFutureEnergy(f, g, h, eta=eta, degree=d);
    
    % Initial condition from Mark Embree's talk
    L = 30; x = linspace(0, L, numEl + 1).';
    initialCondition = x0 * x .* (x - L) .* (x - L / 2);
    
    initialCondition = initialCondition(2:end - 1);
    
    wzInit = 0.5 * kronPolyEval(w, initialCondition, d);
    fprintf('%12.6e    \n', wzInit);
    energies = [energies, wzInit];
end

if exportData
    fprintf('Writing data to plots/example8_convergenceDataWRTd.dat \n')
    fileID = fopen('plots/example8_convergenceDataWRTd.dat', 'w');
    fprintf(fileID, '# Table II Data\n');
    fprintf(fileID, '# finite element heat equation model, convergence and scalability results; n=%d \n', numEl-1);
    %print the header
    fprintf(fileID, 'd & E_d^+(x_0)      \n');
    for i = 1:length(numEls)
        fprintf(fileID, '%5d       &%5d & %10.4e    & %8.2e  & %12.6e   &  %2.2f  &  %2.2f \n', numEls(i), numEls(i) - 1, nd(i), times(i), energies(i), a, d);
    end
    fclose(fileID);
end


end
