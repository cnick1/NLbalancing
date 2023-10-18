%  example KS, feedback control of the discretized periodic
%  Kuramoto-Sivashinsky equations using distributed control.
%
%  These values match those in Section IV.D of the reference
%
%   Reference: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki, “Nonlinear
%               balanced truncation: Part 1—computing energy functions,” arXiv,
%               Dec. 2022. doi: 10.48550/ARXIV.2209.07645
%              [2] B. Kramer, S. Gugercin, and J. Borggaard, “Nonlinear balanced
%               truncation: Part 2—model reduction on manifolds,” arXiv, Feb. 2023.
%               doi: 10.48550/ARXIV.2302.02036
%
%  Part of the NLbalancing and QQR repositories.
%%

eta = 0.1;

m = 5;
p = 2;
L = 13.0291;

z_factor = 0.1; % scale factor on the initial conditions.

%%
%  Computational performance of the energy function approximations.
%  Since the initial times are so short, we average nTest times
%
%  This builds TABLE V
%
fprintf('Table V Data\n')
nTest = 10;
degree = 3;

for n = [8, 16, 32] % number of elements, the actual number of states is 2*n
    fprintf('%d & ', 2 * n)
    [f, g, h, zInit] = getSystem4(n, m, p, 1 / L ^ 2);
    zInit = z_factor * zInit;

    tic; for i = 1:nTest, [w] = approxFutureEnergyQB(f, f{2}, g, h, eta, degree); end, tt = toc / nTest;
    fprintf('%10.4e & ', length(w{degree}))
    fprintf('%8.2e & ', tt)

    for d = 2:degree, w{d} = w{d}.'; end
    wzInit = 0.5 * kronPolyEval(w, zInit, degree);
    fprintf('%12.6e \\\\ \n', wzInit)
end

%%
nTest = 1;
degree = 3;

for n = [64, 128, 256, 512]
    fprintf('%d & ', 2 * n)
    [f, g, h, zInit] = getSystem4(n, m, p, 1 / L ^ 2);
    zInit = z_factor * zInit;

    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;
    fprintf('%10.4e & ', length(w{degree}))
    fprintf('%8.2e & ', tt)

    for d = 2:degree, w{d} = w{d}.'; end
    wzInit = 0.5 * kronPolyEval(w, zInit, degree);
    fprintf('%12.6e \\\\ \n', wzInit)
end

fprintf('\n\n')
%%
%  Test convergence of the energy function at a "point" with increasing
%  degree...
%
%  This builds TABLE VI
%
fprintf('Table VI Data\n')

n = 8; % actual system size is double this
[f, g, h, zInit] = getSystem4(n, m, p, 1 / L ^ 2);
zInit = z_factor * zInit;
for degree = [2, 3, 4, 5, 6] % only the last one is needed, but use this for timing
    fprintf('%d & ', degree)

    tic;
    for kount = 1:2 ^ (8 - degree)
        [v] = approxPastEnergy(f, g, h, eta, degree, false);
    end
    tt = toc / 2 ^ (8 - degree);
    for d = 2:degree, v{d} = v{d}.'; end
    vzInit = 0.5 * kronPolyEval(v, zInit, degree);
    fprintf('%13.7e ', vzInit)
    fprintf('(%8.2e) & ', tt)

    tic;
    for kount = 1:2 ^ (8 - degree)
        [w] = approxFutureEnergy(f, g, h, eta, degree, false);
    end
    tt = toc / 2 ^ (8 - degree);
    for d = 2:degree, w{d} = w{d}.'; end
    wzInit = 0.5 * kronPolyEval(w, zInit, degree);
    fprintf('%13.7e ', wzInit)
    fprintf('(%8.2e) \\\\ \n', tt)
end
save('KSenergyFcns.mat')
