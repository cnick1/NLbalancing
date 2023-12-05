opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [4, Inf];
opts.Delimiter = "&";
opts.VariableNames = ["ns", "energyTimings", "transformationTimings"];
opts.VariableTypes = ["double", "double", "double"];

% Import the data
tbl = readtable("plots\example17_timingData1.dat", opts);

% Calculate the average of entries corresponding to each value of n
n = unique(tbl.ns); % Vector of n values
avg_energy_function_cpu_sec = zeros(size(n));
avg_transformation_cpu_sec = zeros(size(n));
for i = 1:length(n)
    % Find indices of entries with same value of n
    idx = find(tbl.ns == n(i));
    
    % Calculate the average of entries corresponding to each value of n
    avg_energy_function_cpu_sec(i) = mean(tbl.energyTimings(idx));
    avg_transformation_cpu_sec(i) = mean(tbl.transformationTimings(idx));
end

if true %exportData
    fprintf('Writing data to plots/example17_averagedTimings.dat \n')
    fileID = fopen('plots/example17_averagedTimings.dat', 'w');
    fprintf(fileID, '# Table I Data\n');
    fprintf(fileID, '# mass spring damper model, timing results; \n');
    %print the header
    fprintf(fileID, '    n   & Energy function CPU-sec   & Transformation CPU-sec \n');
    for i = 1:length(n)
        fprintf(fileID, '%5d   &        %10.4e         &      %10.4e \n', n(i), avg_energy_function_cpu_sec(i), avg_transformation_cpu_sec(i));
    end
    fclose(fileID);
end

loglog(n, avg_energy_function_cpu_sec);
hold on
loglog(n, avg_transformation_cpu_sec);

