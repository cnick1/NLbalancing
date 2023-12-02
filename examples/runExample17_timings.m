function runExample17_timings()
%runExample17_timings Runs the example to test diagonalization.
%
%   Usage:  [v,w] = runExample17_timings()
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 17\n')
eta = 0;


degree = 4;

energyTimings = [];
transformationTimings = [];

ns = [4,8,16,32,64,128];
for n=ns
    [f, g, h] = getSystem17(degree - 1, n / 2);
    
    %  Compute the energy functions
    tic
    [v] = approxPastEnergy(f, g, h, eta, degree, false);
    [w] = approxFutureEnergy(f, g, h, eta, degree, false);
    energyTimings = [energyTimings, toc];
    fprintf("Computing the energy functions took %f seconds. \n", energyTimings(end))
    
    
    %% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
    tic
    [sigmaSquared, Tod] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
    transformationTimings = [transformationTimings, toc];
    fprintf("Input-normal/output-diagonal transformation took %f seconds. \n", transformationTimings(end))
end

loglog(ns,energyTimings)
hold on
loglog(ns,transformationTimings)

end
