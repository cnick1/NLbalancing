function [ns,energyTimings,transformationTimings] = runExample17_timings(degree)
%runExample17_timings Runs the example to test diagonalization.
%
%   Usage:  [v,w] = runExample17_timings()
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 17\n')

if nargin < 1
    degree = 3;
end

if true %exportData
    fprintf('Writing data to plots/example17_timingData_d%d.dat \n', degree)
    fileID = fopen(sprintf('plots/example17_timingData_d%d.dat', degree), 'w');
    fprintf(fileID, '# Table I Data\n');
    fprintf(fileID, '# mass spring damper model, timing results; \n');
    %print the header
    fprintf(fileID, '    n   & Energy function CPU-sec   & Transformation CPU-sec \n');
end

eta = 0;


energyTimings = [];
transformationTimings = [];
i = 0;

if degree == 6
    ns = [4,8,16];
    % ns = [4*ones(1,10),8*ones(1,10),16*ones(1,3)];
elseif degree == 5
    ns = [4,8,16,32];
    % ns = [4*ones(1,10),8*ones(1,10),16*ones(1,3),32];
elseif degree == 4
    ns = [4,8,16,32,64,128];
    % ns = [4*ones(1,10),8*ones(1,10),16*ones(1,10),32*ones(1,3),64*ones(1,3),128];
elseif degree == 3
    ns = [4,8,16,32,64,128, 256, 512];
    % ns = [4*ones(1,10),8*ones(1,10),16*ones(1,10),32*ones(1,3),64*ones(1,3),128, 256, 512];
end

for n=ns
    [f, g, h] = getSystem17(degree - 1, n / 2);
    
    %  Compute the energy functions
    fprintf("Computing energy functions for n=%i, d=%i ... ", n, degree); tic
    [v] = approxPastEnergy(f, g, h, eta, degree, false);
    [w] = approxFutureEnergy(f, g, h, eta, degree, false);
    energyTimings = [energyTimings, toc];
    fprintf("completed in %f seconds. \n", energyTimings(end))
    
    clear f g h % save some ram
    %% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
    tic
    [sigmaSquared, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
    transformationTimings = [transformationTimings, toc];
    fprintf("\n    -> Transformation computed in %f seconds. \n\n", transformationTimings(end))
    
    if true %exportData
        i = i+1;
        fprintf(fileID, '%5d   &        %10.4e         &      %10.4e \n', ns(i), energyTimings(i), transformationTimings(i));
    end
end


if true %exportData
    fclose(fileID);
end

% Calculate the average of entries corresponding to each value of n
n = unique(ns); % Vector of n values
avg_energy_function_cpu_sec = zeros(size(n));
avg_transformation_cpu_sec = zeros(size(n));
for i = 1:length(n)
    % Find indices of entries with same value of n
    idx = find(ns == n(i));
    
    % Calculate the average of entries corresponding to each value of n
    avg_energy_function_cpu_sec(i) = mean(energyTimings(idx));
    avg_transformation_cpu_sec(i) = mean(transformationTimings(idx));
end

if true %exportData
    fprintf('Writing data to plots/example17_averagedTimings_d%d.dat \n', degree)
    fileID = fopen(sprintf('plots/example17_averagedTimings_d%d.dat', degree), 'w');
    fprintf(fileID, '# Table I Data\n');
    fprintf(fileID, '# mass spring damper model, timing results; \n');
    %print the header
    fprintf(fileID, '    n   & Energy function CPU-sec   & Transformation CPU-sec \n');
    for i = 1:length(n)
        fprintf(fileID, '%5d   &        %10.4e         &      %10.4e \n', n(i), avg_energy_function_cpu_sec(i), avg_transformation_cpu_sec(i));
    end
    fclose(fileID);
end

end
